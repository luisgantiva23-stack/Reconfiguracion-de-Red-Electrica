% MÉTODO DE APROXIMACIONES SUCESIVAS
clc
clear
%% Datos del sistema
%        i j Rij(Ohm) Xij(Ohm)
%        i   j     R(Ohm)   X(Ohm)
Lines = [
1   2   0.0922   0.0470;
2   3   0.3930   0.2512;
3   4   0.2661   0.1864;
%4   5   0.1725   0.1435;
5   6   0.1899   0.4114;
1   7   0.2565   0.3624;
7   8   0.1416   0.1359;
%7   9   0.1598   0.1624;
%8   10  0.3158   0.2886;
2   11  0.2689   0.3154;
11  12  0.2416   0.2548;
%11  13  0.3050   0.2159;
12  14  0.3536   0.3333;
%4   15  0.2125   0.3895;

%15  16  0.3654   0.3456;
%5   17  0.2224   0.2550;
%17  18  0.3587   0.3356;
%6   19  0.3025   0.3486;
19  20  0.2986   0.1457;
%19  21  0.1954   0.1854;
%6   22  0.1489   0.1875;
%22  23  0.3589   0.1752;
%22  24  0.2222   0.2758;
2   9   0.1645   0.1325;
%3   9   0.1196   0.1745;
3   13  0.2565   0.3528;
5   20  0.1654   0.1212;
5   22  0.2958   0.2566;

9   10  0.1415   0.1111;
9   15  0.1784   0.2256;
10  16  0.1212   0.1456;
%13  14  0.3021   0.3022;
13  17  0.1229   0.1565;
14  18  0.2560   0.2222;
15  20  0.1475   0.1865;
%16  20  0.2898   0.3022;
17  23  0.1652   0.1895;
%18  23  0.3012   0.3526;
%19  24  0.1238   0.1985;
20  21  0.2547   0.2259;
21  24  0.1331   0.1441
];
%        i Pi(kW) Qi(kvar)
%        i   Pd(kW)   Qd(kvar)
Nodes = [
1   0    0;
2   145  100;
3   175  120;
4   400  325;
5   700  600;
6   325  200;
7   0    -900;
8   500  400;
9   800  750;
10  900  700;
11  850  600;
12  0    -1250;
13  350  300;
14  700  625;
15  600  500;
16  650  600;
17  500  400;
18  600  500;
19  750  500;
20  450  300;
21  0    -1750;
22  600  750;
23  850  800;
24  650  400
];

Vb = 11400; % Voltaje base
Sb = 1000;  % Potencia base kva
Zb = (Vb^2)/(Sb*1000);
% Pasar a pu los datos del sistema
Lines(:,3:4) = Lines(:,3:4)/Zb;
Nodes(:,2:3) = Nodes(:,2:3)/Sb;
N = size(Nodes,1); % Número de nodos
L = size(Lines,1); % Número de líneas
%% Construir la matriz de admitancia nodal (Ybus)
Ybus = zeros(N,N);
for l = 1:L
    k = Lines(l,1); m = Lines(l,2); Zkm = Lines(l,3) + 1i*Lines(l,4);
    Ybus(k,k) = Ybus(k,k) + 1/Zkm;
    Ybus(m,m) = Ybus(m,m) + 1/Zkm;
    Ybus(k,m) = -1/Zkm;
    Ybus(m,k) = -1/Zkm;
end
% Extraer información para el flujo (Se asume slack en el nodo 1)
Yss = Ybus(1,1);
Ysd = Ybus(1,2:N);
Yds = Ybus(2:N,1);
Ydd = Ybus(2:N,2:N);
Zdd = inv(Ydd);
Sd = Nodes(2:N,2) + 1i*Nodes(2:N,3); % P + jQ
Vs = 1; % Vslack = 1<0
Vd0 = ones(N-1,1);
error = 1e-4;
tmax = 100;
% Ciclo para el cálculo del flujo
for t = 1:tmax
    Vd = -Zdd*(Yds*Vs + inv(diag(conj(Vd0)))*conj(Sd)); %#ok<MINV>
    err = max(abs(abs(Vd) - abs(Vd0)));
    if err <= error
        V = [Vs;Vd];
        break
    else
        Vd0 = Vd;
    end
end
% Imprimir resultados
fprintf('Nodo \t Mag. V \t Ang. V\n')
for k = 1:N
    fprintf('%d\t %f\t %f\t\n',k,abs(V(k)),angle(V(k))*180/pi);
end
% Pérdidas totales
Sgt = Vs*conj(Yss*Vs + Ysd*V(2:N,1));
Sdt = sum(Sd);
Sloss = Sgt - Sdt;
% Cálculos para las líneas
fprintf('Envío Recibo\tMag. I \t Ang. I \t Perd (kw)\t Perd. (kvar)\n')
pl = 0;
for l = 1:L
    k = Lines(l,1); m = Lines(l,2); Zkm = Lines(l,3) + 1i*Lines(l,4);
    Il = (V(k)-V(m))/Zkm;
    Pl = (abs(Il))^2*real(Zkm)*Sb;
    pl = Pl + pl;
    Ql = (abs(Il))^2*imag(Zkm)*Sb;
    fprintf('%d\t %d\t \t%f\t %f\t %f\t %f\t\n',k,m,abs(Il),angle(Il)*180/pi,Pl,Ql);
end

Ploss = real(transpose(V)*conj(Ybus*V)*Sb);
fprintf('Perdidas totales: %f\t\n', Ploss);
%% Factor de potencia de la subestacion
Pgt = real(Vs*conj(Yss*Vs + Ysd*V(2:N,1)));
Fpsub = Pgt/abs(Sgt);

%% Regulacion de tension
for k = 2:N
    RV_k = (abs(Vs) - abs(V(k))) / abs(Vs) * 100;
    if abs(V(k)) >= 0.95 && abs(V(k)) <= 1.05
        estado = 'OK';
    else
        estado = 'FUERA';
    end
    fprintf('%d\t\t %.4f\t\t %.2f%%\t\t %s\n', k, abs(V(k)), RV_k, estado);
end
% Regulacion Mas alta
Vmag = abs(V);
[Vmin, nodo_Vmin] = min(Vmag(2:N));   % Excluye el slack
nodo_Vmin = nodo_Vmin + 1;            % Índice real en el sistema
[Vmax_carga, nodo_Vmax] = max(Vmag(2:N));
nodo_Vmax = nodo_Vmax + 1;
RV = (abs(Vs) - Vmin) / abs(Vs) * 100; % Regulación de tensión (%)
fprintf('            RESUMEN DEL CASO BASE\n')
fprintf('====================================================\n')
fprintf('Pérdidas totales activas:       %.4f kW\n', Ploss)
fprintf('Regulación de tensión maxima:   %.2f %%\n', RV)
fprintf('Factor de potencia subestación: %.4f\t\n', Fpsub);
fprintf('====================================================\n')
