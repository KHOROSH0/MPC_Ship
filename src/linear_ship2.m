clear all;
close all;

% шаг дискретизации
T = 1;
T0 = 0;
T_final = 100

% начальные условия
x0 = [0.01; 0.01; 0; 0.03]*0;
r = 20*pi/180;

%ограничения на управление
u_d = 3*pi/180;
delta_d = 35*pi/180;

% время движения
t = T0:T:T_final;

% горизонт прогноза
p = 50;

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
B = [0;0;0;1]
C = [0,0,1,0];
C_r = [0,0,0,1];
%C = eye(4,4);
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
[K,S,aaaa] = dlqr(A_d,B_d,Q,R);
K = -K;
eig(A_d+B_d*K);

% моделирование LQR-регулятора

A_cl = (A_d+B_d*K);
B_cl = zeros(4,1);
C_cl = C;

sys_cl = ss(A_cl,B_cl,C_cl,D_d,T);


% график модели замкнутой LQR-регулятором 
%initial(sys_cl, x0);
%grid on;


% MPC подход с разными p (горизонт прогноза)

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
K_mpc = -inv((M1')*Q_mpc*M1 + R_mpc)*(M1')*Q_mpc*L1;
T_mpc = inv((M1')*Q_mpc*M1 + R_mpc)*(M1')*Q_mpc;

%K_mpc = K_mpc(1,:)
%T_mpc = T_mpc(1,:)


% Второе задание 

% Составим допустимое множество для регулятора
x_ = x0;
r_ = r;
by = [delta_d];
bu = [u_d];
for i=2:p
    r_ = [r_;r];
    by = [by; delta_d];
    bu = [bu;u_d];
end
by = [by;by];
bu = [bu;bu];
% Ограничение на величину управления
Au = [eye(p,p); -eye(p,p)];

% Ограничение на контролируемую переменную в нашем случае дельта
%построим матрицы M1 и L1 для ограничения 

M1_d = C_r*B_d;
M_str_d = B_d;
L1_d = C_r*A_d;
for i =2:p
    L1_d = [L1_d;C_r*(A_d^i)];
    M_str_d = [A_d*M_str_d,B_d];
    [rows, cols] = size(M1_d);
    M1_d = [M1_d, zeros(rows,1); C_r*M_str_d];
end
Ay = [M1_d; -M1_d];
bx = [-L1_d; L1_d];


% составим общую систему ограничений
A_big = [Au; Ay];
b = [bu;by] + [zeros(size(bu)); bx*x_];


% моделирование с вычислением управления на каждом шаге с учетом
% ограничений
% первый прогон с нулевыми управлениями, надеюсь MPC дальше исправится
% составим вектор 
g = L1*x_ - r_
H = (M1'*Q_mpc*M1 + R_mpc)
H = (H + H')/2;
f = M1'*Q_mpc*g;
%H = (H + H')/2
J = @(u) (u'*H*u) + 2*f'*u + g'*Q_mpc*g;
u_ = [0];
x_array3 = [];
x_array4 = [];
u_array = [];
opt = optimset('Display', 'iter', 'TolFun', 1e-15, 'TolX', 1e-15, 'MaxIter', 100);

for i=T0:T:T_final
    [u_, x, a] = quadprog(H,f,A_big,b,[],[],[],[],[],opt);
    value = (A_big*u_ - b) < -1e-10;
    if sum(value) == 4*p
        u_ = K_mpc*x_ + T_mpc*r_;
    end
    x_ = A_d*x_ + B_d*u_(1);
    b = [bu;by] + [zeros(size(bu)); bx*x_];
    g = L1*x_ - r_;
    f = M1'*Q_mpc*g;
    x_array3 =[x_array3, x_(3)];
    x_array4 =[x_array4, x_(4)];
    u_array = [u_array, u_(1)];
end
subplot(1,3,1)
plot(t, x_array3*180/pi)
title('angle of ship')
xlabel('time')
ylabel('\phi')
grid on
subplot(1,3,2)
plot(t, x_array4*180/pi)
title('angle of wheel ')
xlabel('time')
ylabel('\delta')
grid on
subplot(1,3,3)
plot(t,u_array*180/pi)
title('velocity of wheel (control)')
xlabel('time')
ylabel('u')
grid on
