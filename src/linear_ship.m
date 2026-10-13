clear all;
close all;

% шаг дискретизации
T = 1;
% начальные условия

x0 = [0.1, 0.2, 0, 0.3];

% время движения
t = 0:T:5;

% Определение переменных (по смыслу из текста)
% beta (β) - угол дрейфа
% omega (ω) - угловая скорость по курсу
% phi (φ) - угол рыскания
% delta (δ) - угол поворота рулей

% Задание коэффициентов
a11 = -0.0454;
a12 = 0.560;
a21 = 0.0267;
a22 = -0.408;

b1 = -0.0132;
b2 = -0.00742;

h1 = -0.0648;
h2 = -0.00456;

% Составление матриц

A = [a11, a12, 0, b1;
    a21, a22, 0, b2;
    0, 1, 0, 0;
    0, 0, 0, 0];
B1 = [0;0;0;1]
B = B1
C = [0,0,1,0];
D = 0

sys = ss(A,B,C,D)


% Дискретная система
%A_d = A*T + eye(4,4);
%B_d = T*B;
%C_d = C;
%D_d = zeros(1,1);
%ship_sys = ss(A_d, B_d,C_d,D_d,T);
% дисретизация C2D

sys_d = c2d(sys,T)

A_d = sys_d.A;
B_d = sys_d.B;
C_d = sys_d.C;
D_d = sys_d.D;


% дискретный функционал
R = eye(1,1)
Q = [0, 0, 0, 0; 
     0, 0, 0, 0; 
     0, 0, 1, 0; 
     0, 0, 0, 0];

% LQR-регулятор
[K,S,eye] = dlqr(A_d,B_d,Q,R);
K = -K;
eig(A_d+B_d*K)

% моделирование LQR-регулятора

A_cl = (A_d+B_d*K);
B_cl = zeros(4,1);
C_cl = C;

sys_cl = ss(A_cl,B_cl,C_cl,D_d,T);


% график модели замкнутой LQR-регулятором 
%initial(sys_cl, x0);
%grid on;


% MPC подход с разными p (горизонт прогноза)
p = 50;

% алгоритм построения
% строим матрицу L1 и M1
M1 = C*B_d;
M_str = B_d;
L1 = C*A_d;
for i =2:p
    L1 = [L1;C*(A_d^i)];
    M_str = [A_d*M_str,B_d];
    [rows, cols] = size(M1);
    M1 = [M1, zeros(rows,1); C*M_str];
end


% вычисляем K и T
% строим Q_mpc и R_mpc
Q_mpc = [1];
Q = [1];
R_mpc = R;

[rowsQ, colsQ] = size(Q);
[rowsR, colsR] = size(R);

for i =2:p

    [rowsQ_mpc, colsQ_mpc] = size(Q_mpc);
    [rowsR_mpc, colsR_mpc] = size(R_mpc);
    Q_mpc = [Q_mpc,zeros(rowsQ_mpc,colsQ);zeros(rowsQ,colsQ_mpc),Q];
    R_mpc = [R_mpc,zeros(rowsR_mpc,colsR);zeros(rowsR,colsR_mpc),R];
end

% строим K_ и T_
K_mpc = -inv((M1.')*Q_mpc*M1 + R_mpc)*(M1.')*Q_mpc*L1;
T_mpc = inv((M1.')*Q_mpc*M1 + R_mpc)*(M1.')*Q_mpc;



%sum = B_d
%for i = 2:p
%   sum = [A_d*sum, B_d];
%end
%result = sum*K_mpc;
disp("K_mpc")
disp(eig(A_d + B_d*K_mpc(1,:)))
disp(eig(A_d+B_d*K))
disp(K_mpc(1,:))
disp(K)


    






