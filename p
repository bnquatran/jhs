%%
clear all
clc

%%
global rho cp T_air h k T_b
syms x spacing N lambda1 lambda2 lambda3 lambda4 s1 s2 s3

% units: W, J, m, °C, kg
rho = 1.225;
cp = 1000;
T_air = 25;
h = 100;
k = 287;
T_b = 85;
L_cpu = 70 * 10^-3;
q_cpu_dot = 100;

P = 4 * x;
A_cpu = L_cpu^2;
m = sqrt((h * P) / (k * A_cpu));
A_b = (N * x + (N - 1) * spacing)^2 - N * A_cpu;
theta_b = T_b - T_air;
M = theta_b * sqrt(h * P * k * A_cpu);
q_fin_dot = M * (sinh(m * L_cpu) + h / (m * k) * cosh(m * L_cpu)) / (cosh(m * L_cpu) + h / (m * k) * sinh(m * L_cpu));

N = solve(q_cpu_dot == N * q_fin_dot + h * A_b * theta_b, N);

phi1 = ((L_cpu + spacing) / (x + spacing))^2;


L = N - lambda1 * phi1;

e1 = diff(L, x) == 0;
e2 = diff(L, spacing) == 0;
e3 = diff(L, lambda1) == 0;

sol = solve([e1, e2, e3], [x, spacing, lambda1]);

%% 

sol.x
sol.spacing
sol.lambda1
