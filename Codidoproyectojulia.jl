using JuMP
using Ipopt
using Juniper
using HiGHS
using LinearAlgebra
using Printf

function es_radial(lines_idx::Vector{Int}, N::Int, from::Vector{Int}, to::Vector{Int}, slack::Int)
    if length(lines_idx) != N - 1
        return false
    end

    adj = [Int[] for _ in 1:N]
    for l in lines_idx
        i = from[l]
        j = to[l]
        push!(adj[i], j)
        push!(adj[j], i)
    end

    visitado = falses(N)
    padre = zeros(Int, N)

    cola = [slack]
    visitado[slack] = true

    while !isempty(cola)
        nodo = popfirst!(cola)
        for vecino in adj[nodo]
            if !visitado[vecino]
                visitado[vecino] = true
                padre[vecino] = nodo
                push!(cola, vecino)
            elseif vecino != padre[nodo]
                return false
            end
        end
    end

    return all(visitado)
end

function flujo_potencia_sa(
    lines_idx::Vector{Int},
    N::Int,
    branch_data_all,
    Zb,
    from::Vector{Int},
    to::Vector{Int},
    R::Vector{Float64},
    X::Vector{Float64},
    Sd::Vector{ComplexF64},
    Sb_kVA::Float64,
    slack::Int;
    tol=1e-6,
    maxiter=200
)
    lines_pu = hcat(
        branch_data_all[lines_idx, 1:2],
        branch_data_all[lines_idx, 3:4] ./ Zb
    )

    Ybus = zeros(ComplexF64, N, N)

    for l in 1:size(lines_pu, 1)
        i = Int(lines_pu[l, 1])
        j = Int(lines_pu[l, 2])
        z = lines_pu[l, 3] + 1im * lines_pu[l, 4]
        y = 1 / z

        Ybus[i, i] += y
        Ybus[j, j] += y
        Ybus[i, j] -= y
        Ybus[j, i] -= y
    end

    Yds = Ybus[2:N, 1]
    Ydd = Ybus[2:N, 2:N]

    Zdd = try
        inv(Ydd)
    catch
        return (ok=false,)
    end

    Vs = 1.0 + 0im
    Vd = ones(ComplexF64, N - 1)

    convergio = false
    for _ in 1:maxiter
        Dinv = Diagonal(1.0 ./ conj.(Vd))
        Vd_new = -Zdd * (Yds * Vs + Dinv * conj.(Sd[2:N]))

        if maximum(abs.(abs.(Vd_new) .- abs.(Vd))) < tol
            Vd = Vd_new
            convergio = true
            break
        end
        Vd = Vd_new
    end

    if !convergio
        return (ok=false,)
    end

    V = [Vs; Vd]

    Ploss = 0.0
    Qloss = 0.0

    for l in lines_idx
        i = from[l]
        j = to[l]
        z = R[l] + 1im * X[l]
        Iij = (V[i] - V[j]) / z

        Ploss += real(z) * abs2(Iij) * Sb_kVA
        Qloss += imag(z) * abs2(Iij) * Sb_kVA
    end

    Isub = sum(Ybus[slack, m] * V[m] for m in 1:N)
    Ssub = V[slack] * conj(Isub) * Sb_kVA

    Vmag = abs.(V)
    Vmin = minimum(Vmag)
    Vmax = maximum(Vmag)
    nodo_min = argmin(Vmag)

    return (
        ok=true,
        V=V,
        Vmag=Vmag,
        Ploss=Ploss,
        Qloss=Qloss,
        Pg=real(Ssub),
        Qg=imag(Ssub),
        Ssub=abs(Ssub),
        fp=abs(Ssub) > 1e-9 ? real(Ssub) / abs(Ssub) : 1.0,
        Vmin=Vmin,
        Vmax=Vmax,
        nodo_min=nodo_min
    )
end

function orientar_arbol(
    lines_idx::Vector{Int},
    N::Int,
    from::Vector{Int},
    to::Vector{Int},
    slack::Int,
    L::Int
)
    adj = [Int[] for _ in 1:N]
    adj_line = [Int[] for _ in 1:N]

    for l in lines_idx
        i = from[l]
        j = to[l]
        push!(adj[i], j)
        push!(adj[j], i)
        push!(adj_line[i], l)
        push!(adj_line[j], l)
    end

    visitado = falses(N)
    profundidad = zeros(Float64, N)
    yF0 = zeros(Float64, L)
    yR0 = zeros(Float64, L)

    cola = [slack]
    visitado[slack] = true

    while !isempty(cola)
        u = popfirst!(cola)
        for k in eachindex(adj[u])
            v = adj[u][k]
            l = adj_line[u][k]

            if !visitado[v]
                visitado[v] = true
                profundidad[v] = profundidad[u] + 1
                if from[l] == u && to[l] == v
                    yF0[l] = 1.0
                else
                    yR0[l] = 1.0
                end
                push!(cola, v)
            end
        end
    end

    return yF0, yR0, profundidad
end

function heuristica_warmstart(
    lines_base::Vector{Int},
    N::Int,
    branch_data_all,
    Zb,
    from::Vector{Int},
    to::Vector{Int},
    R::Vector{Float64},
    X::Vector{Float64},
    Sd::Vector{ComplexF64},
    Sb_kVA::Float64,
    slack::Int
)
    lines_actual = copy(lines_base)

    eval_actual = flujo_potencia_sa(
        lines_actual, N, branch_data_all, Zb, from, to, R, X, Sd, Sb_kVA, slack
    )

    if !get(eval_actual, :ok, false)
        error("No se pudo evaluar el caso base con aproximación sucesiva.")
    end

    mejoro = true

    while mejoro
        mejoro = false
        mejor_eval = eval_actual
        mejor_config = copy(lines_actual)

        for l_entra in 24:41
            if l_entra in lines_actual
                continue
            end

            for l_sale in lines_actual
                config_prueba = sort(vcat(filter(x -> x != l_sale, lines_actual), [l_entra]))

                if !es_radial(config_prueba, N, from, to, slack)
                    continue
                end

                eval_prueba = flujo_potencia_sa(
                    config_prueba, N, branch_data_all, Zb, from, to, R, X, Sd, Sb_kVA, slack
                )

                if !get(eval_prueba, :ok, false)
                    continue
                end

                if eval_prueba.Ploss < mejor_eval.Ploss - 1e-6
                    mejor_eval = eval_prueba
                    mejor_config = copy(config_prueba)
                    mejoro = true
                end
            end
        end

        if mejoro
            lines_actual = mejor_config
            eval_actual = mejor_eval
        end
    end

    return lines_actual, eval_actual
end

function resolver_minlp(
    warm_lines::Vector{Int},
    N::Int,
    L::Int,
    from::Vector{Int},
    to::Vector{Int},
    branch_data_all,
    R::Vector{Float64},
    X::Vector{Float64},
    Pd::Vector{Float64},
    Qd::Vector{Float64},
    Sd::Vector{ComplexF64},
    Zb::Float64,
    Sb_kVA::Float64,
    slack::Int,
    VMIN2::Float64,
    VMAX2::Float64,
    Pmax::Float64,
    Qmax::Float64,
    I2max::Float64,
    Mbig::Float64;
    max_changes::Union{Nothing, Int}=nothing
)
    warm_eval = flujo_potencia_sa(
        warm_lines, N, branch_data_all, Zb, from, to, R, X, Sd, Sb_kVA, slack
    )
    if !get(warm_eval, :ok, false)
        error("No se pudo evaluar el warm start.")
    end

    yF0, yR0, depth0 = orientar_arbol(warm_lines, N, from, to, slack, L)
    z0 = zeros(Float64, L)
    z0[warm_lines] .= 1.0

    nl_solver = optimizer_with_attributes(
        Ipopt.Optimizer,
        "print_level" => 0,
        "max_iter" => 2500,
        "tol" => 1e-6,
        "acceptable_tol" => 1e-4,
        "acceptable_iter" => 10
    )

    mip_solver = optimizer_with_attributes(
        HiGHS.Optimizer,
        "output_flag" => false
    )

    model = Model(
        optimizer_with_attributes(
            Juniper.Optimizer,
            "nl_solver" => nl_solver,
            "mip_solver" => mip_solver,
            "time_limit" => 240.0
        )
    )

    set_silent(model)

    @variable(model, z[1:L], Bin)
    @variable(model, yF[1:L], Bin)
    @variable(model, yR[1:L], Bin)

    @variable(model, 0 <= u[1:N] <= N - 1)
    @variable(model, VMIN2 <= V2[1:N] <= VMAX2)

    @variable(model, 0 <= PF[1:L] <= Pmax)
    @variable(model, 0 <= PR[1:L] <= Pmax)
    @variable(model, -Qmax <= QF[1:L] <= Qmax)
    @variable(model, -Qmax <= QR[1:L] <= Qmax)
    @variable(model, 0 <= I2F[1:L] <= I2max)
    @variable(model, 0 <= I2R[1:L] <= I2max)

    @variable(model, 0 <= Pg0 <= Pmax + 2.0)
    @variable(model, -Qmax - 2.0 <= Qg0 <= Qmax + 2.0)

    fix(V2[slack], 1.0; force=true)
    fix(u[slack], 0.0; force=true)

    for l in 1:L
        set_start_value(z[l], z0[l])
        set_start_value(yF[l], yF0[l])
        set_start_value(yR[l], yR0[l])
        set_start_value(PF[l], 0.0)
        set_start_value(PR[l], 0.0)
        set_start_value(QF[l], 0.0)
        set_start_value(QR[l], 0.0)
        set_start_value(I2F[l], 0.0)
        set_start_value(I2R[l], 0.0)
    end

    for n in 1:N
        set_start_value(u[n], depth0[n])
        set_start_value(V2[n], abs(warm_eval.V[n])^2)
    end
    set_start_value(Pg0, warm_eval.Pg / Sb_kVA)
    set_start_value(Qg0, warm_eval.Qg / Sb_kVA)

    @constraint(model, [l=1:L], yF[l] + yR[l] == z[l])
    @constraint(model, sum(z[l] for l in 1:L) == N - 1)

    for n in 1:N
        expr_in = AffExpr()
        for l in 1:L
            if to[l] == n
                add_to_expression!(expr_in, yF[l])
            end
            if from[l] == n
                add_to_expression!(expr_in, yR[l])
            end
        end
        if n == slack
            @constraint(model, expr_in == 0)
        else
            @constraint(model, expr_in == 1)
        end
    end

    for l in 1:L
        i = from[l]
        j = to[l]
        @constraint(model, u[j] >= u[i] + 1 - N * (1 - yF[l]))
        @constraint(model, u[i] >= u[j] + 1 - N * (1 - yR[l]))
    end

    if max_changes !== nothing
        expr_change = AffExpr()
        for l in 1:L
            if z0[l] > 0.5
                add_to_expression!(expr_change, 1.0)
                add_to_expression!(expr_change, -1.0, z[l])
            else
                add_to_expression!(expr_change, z[l])
            end
        end
        @constraint(model, expr_change <= max_changes)
    end

    @constraint(model, [l=1:L], PF[l] <= Pmax * yF[l])
    @constraint(model, [l=1:L], PR[l] <= Pmax * yR[l])

    @constraint(model, [l=1:L], QF[l] <=  Qmax * yF[l])
    @constraint(model, [l=1:L], QF[l] >= -Qmax * yF[l])
    @constraint(model, [l=1:L], QR[l] <=  Qmax * yR[l])
    @constraint(model, [l=1:L], QR[l] >= -Qmax * yR[l])

    @constraint(model, [l=1:L], I2F[l] <= I2max * yF[l])
    @constraint(model, [l=1:L], I2R[l] <= I2max * yR[l])

    for n in 1:N
        exprP = AffExpr()
        exprQ = AffExpr()

        for l in 1:L
            if to[l] == n
                add_to_expression!(exprP, PF[l])
                add_to_expression!(exprP, -R[l], I2F[l])
                add_to_expression!(exprQ, QF[l])
                add_to_expression!(exprQ, -X[l], I2F[l])
            end
            if from[l] == n
                add_to_expression!(exprP, PR[l])
                add_to_expression!(exprP, -R[l], I2R[l])
                add_to_expression!(exprQ, QR[l])
                add_to_expression!(exprQ, -X[l], I2R[l])
            end

            if from[l] == n
                add_to_expression!(exprP, -1.0, PF[l])
                add_to_expression!(exprQ, -1.0, QF[l])
            end
            if to[l] == n
                add_to_expression!(exprP, -1.0, PR[l])
                add_to_expression!(exprQ, -1.0, QR[l])
            end
        end

        if n == slack
            add_to_expression!(exprP, 1.0, Pg0)
            add_to_expression!(exprQ, 1.0, Qg0)
        end

        add_to_expression!(exprP, -Pd[n])
        add_to_expression!(exprQ, -Qd[n])

        @constraint(model, exprP == 0)
        @constraint(model, exprQ == 0)
    end

    for l in 1:L
        i = from[l]
        j = to[l]
        z2 = R[l]^2 + X[l]^2

        @constraint(model,
            V2[j] - V2[i] + 2 * (R[l] * PF[l] + X[l] * QF[l]) - z2 * I2F[l] <= Mbig * (1 - yF[l]))
        @constraint(model,
            V2[j] - V2[i] + 2 * (R[l] * PF[l] + X[l] * QF[l]) - z2 * I2F[l] >= -Mbig * (1 - yF[l]))

        @constraint(model,
            V2[i] - V2[j] + 2 * (R[l] * PR[l] + X[l] * QR[l]) - z2 * I2R[l] <= Mbig * (1 - yR[l]))
        @constraint(model,
            V2[i] - V2[j] + 2 * (R[l] * PR[l] + X[l] * QR[l]) - z2 * I2R[l] >= -Mbig * (1 - yR[l]))
    end

    @NLconstraint(model, [l=1:L], PF[l]^2 + QF[l]^2 <= V2[from[l]] * I2F[l] + 1e-8)
    @NLconstraint(model, [l=1:L], PR[l]^2 + QR[l]^2 <= V2[to[l]]   * I2R[l] + 1e-8)

    @objective(model, Min, sum(R[l] * (I2F[l] + I2R[l]) for l in 1:L))

    optimize!(model)

    status = string(termination_status(model))

    zval = try
        value.(z)
    catch
        error("El MINLP terminó sin una solución legible. Estado: $status")
    end

    selected = findall(v -> v > 0.5, zval)

    if length(selected) != N - 1 || !es_radial(selected, N, from, to, slack)
        return (
            ok=false,
            status=status,
            selected=selected
        )
    end

    final_eval = flujo_potencia_sa(
        selected, N, branch_data_all, Zb, from, to, R, X, Sd, Sb_kVA, slack
    )

    if !get(final_eval, :ok, false)
        return (
            ok=false,
            status=status,
            selected=selected
        )
    end

    return (
        ok=true,
        status=status,
        selected=selected,
        eval=final_eval,
        model_obj=objective_value(model) * Sb_kVA
    )
end

function run_reconfiguracion_minlp()
    Vb_kV = 11.4
    Sb_kVA = 10000.0
    Zb = (Vb_kV * 1000)^2 / (Sb_kVA * 1000)
    slack = 1

    branch_data_all = [
        1   2   0.0922  0.0470
        2   3   0.3930  0.2512
        3   4   0.2661  0.1864
        4   5   0.1725  0.1435
        5   6   0.1899  0.4114
        1   7   0.2565  0.3624
        7   8   0.1416  0.1359
        7   9   0.1598  0.1624
        8   10  0.3158  0.2886
        2   11  0.2689  0.3154
        11  12  0.2416  0.2548
        11  13  0.3050  0.2159
        12  14  0.3536  0.3333
        4   15  0.2125  0.3895
        15  16  0.3654  0.3456
        5   17  0.2224  0.2550
        17  18  0.3587  0.3356
        6   19  0.3025  0.3486
        19  20  0.2986  0.1457
        19  21  0.1954  0.1854
        6   22  0.1489  0.1875
        22  23  0.3589  0.1752
        22  24  0.2222  0.2758
        2   9   0.1645  0.1325
        3   9   0.1196  0.1745
        3   13  0.2565  0.3528
        5   20  0.1654  0.1212
        5   22  0.2958  0.2566
        9   10  0.1415  0.1111
        9   15  0.1784  0.2256
        10  16  0.1212  0.1456
        13  14  0.3021  0.3022
        13  17  0.1229  0.1565
        14  18  0.2560  0.2222
        15  20  0.1475  0.1865
        16  20  0.2898  0.3022
        17  23  0.1652  0.1895
        18  23  0.3012  0.3526
        19  24  0.1238  0.1985
        20  21  0.2547  0.2259
        21  24  0.1331  0.1441
    ]

    node_data = [
        1   0      0
        2   145    100
        3   175    120
        4   400    325
        5   700    600
        6   325    200
        7   0     -900
        8   500    400
        9   800    750
        10  900    700
        11  850    600
        12  0     -1250
        13  350    300
        14  700    625
        15  600    500
        16  650    600
        17  500    400
        18  600    500
        19  750    500
        20  450    300
        21  0     -1750
        22  600    750
        23  850    800
        24  650    400
    ]

    N = size(node_data, 1)
    L = size(branch_data_all, 1)

    from = Int.(branch_data_all[:, 1])
    to   = Int.(branch_data_all[:, 2])

    R = Float64.(branch_data_all[:, 3] ./ Zb)
    X = Float64.(branch_data_all[:, 4] ./ Zb)

    Pd = Float64.(node_data[:, 2] ./ Sb_kVA)
    Qd = Float64.(node_data[:, 3] ./ Sb_kVA)
    Sd = ComplexF64.(Pd .+ 1im .* Qd)

    lines_base = collect(1:23)

    VMIN = 0.90
    VMAX = 1.05
    VMIN2 = VMIN^2
    VMAX2 = VMAX^2

    Ptot = sum(Pd)
    Qtot = sum(abs.(Qd))
    Pmax = Ptot + 1.0
    Qmax = Qtot + 1.0
    I2max = (Pmax^2 + Qmax^2) / VMIN2
    Mbig = (VMAX2 - VMIN2) +
           2 * maximum(abs.(R)) * Pmax +
           2 * maximum(abs.(X)) * Qmax +
           maximum(R.^2 .+ X.^2) * I2max + 1.0

    # Warm start interno, no se imprime
    warm_lines, _ = heuristica_warmstart(lines_base, N, branch_data_all, Zb, from, to, R, X, Sd, Sb_kVA, slack)

    cambios_vecindario = [0, 2, 4]
    candidatos = []

    for K in cambios_vecindario
        sol = resolver_minlp(
            warm_lines, N, L, from, to, branch_data_all, R, X, Pd, Qd, Sd,
            Zb, Sb_kVA, slack, VMIN2, VMAX2, Pmax, Qmax, I2max, Mbig;
            max_changes=K
        )

        if sol.ok
            push!(candidatos, (K=K, status=sol.status, selected=sol.selected, eval=sol.eval))
        end
    end

    if isempty(candidatos)
        error("No hubo soluciones MINLP factibles en los vecindarios evaluados.")
    end

    idx_mejor = argmin([c.eval.Ploss for c in candidatos])
    mejor = candidatos[idx_mejor]

    entran = setdiff(mejor.selected, lines_base)
    salen = setdiff(lines_base, mejor.selected)

    println("\n============================================================")
    println(" SOLUCIÓN FINAL MINLP - SISTEMA DE 24 NODOS")
    println("============================================================")
    @printf("Estado del solver:          %s\n", mejor.status)
    @printf("Vecindario ganador K:       %d\n", mejor.K)
    @printf("Pérdidas activas:           %.3f kW\n", mejor.eval.Ploss)
    @printf("Pérdidas reactivas:         %.3f kvar\n", mejor.eval.Qloss)
    @printf("P subestación:              %.3f kW\n", mejor.eval.Pg)
    @printf("Q subestación:              %.3f kvar\n", mejor.eval.Qg)
    @printf("S subestación:              %.3f kVA\n", mejor.eval.Ssub)
    @printf("Factor de potencia:         %.4f\n", mejor.eval.fp)
    @printf("Tensión mínima:             %.5f pu (Nodo %d)\n", mejor.eval.Vmin, mejor.eval.nodo_min)
    @printf("Tensión máxima:             %.5f pu\n", mejor.eval.Vmax)
    @printf("Líneas activas finales:     %d / %d\n", length(mejor.selected), L)
    @printf("Sistema radial:             %s\n", es_radial(mejor.selected, N, from, to, slack) ? "SÍ" : "NO")

    println("\nCONFIGURACIÓN FINAL DE LÍNEAS ACTIVAS:")
    println("No.   L#    Nodo i -> Nodo j")
    for (k, l) in enumerate(sort(mejor.selected))
        @printf("%2d    L%-2d   %2d -> %2d\n", k, l, from[l], to[l])
    end

    println("\nLÍNEAS QUE ENTRAN A LA CONFIGURACIÓN FINAL:")
    if isempty(entran)
        println("Ninguna")
    else
        println("L#   Nodo i -> Nodo j")
        for l in sort(entran)
            @printf("L%-2d   %2d -> %2d\n", l, from[l], to[l])
        end
    end

    println("\nLÍNEAS QUE SALEN DE LA CONFIGURACIÓN ORIGINAL:")
    if isempty(salen)
        println("Ninguna")
    else
        println("L#   Nodo i -> Nodo j")
        for l in sort(salen)
            @printf("L%-2d   %2d -> %2d\n", l, from[l], to[l])
        end
    end

    println("\nPERFIL DE TENSIÓN FINAL:")
    println("Nodo    |V|(pu)")
    for n in 1:N
        @printf("%2d      %.5f\n", n, mejor.eval.Vmag[n])
    end

    println("============================================================")
end

run_reconfiguracion_minlp()