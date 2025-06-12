clear; clc; close all;

% ASEN 5044
% Midterm 1, 10/10/2024
% Author - Jash Bhalavat

%% Problem 1

% Given constants
g = 9.81;
l = 1.85;
m = 2;
M = 4;

% Load data and assign variables
data = load("midterm1problem1data.mat");
Kc = data.Kc;
thist = data.thist;
yNLhist = data.yNLhist;

% Part c
delta_t = 0.05;

% CT LTI matrices as computed in part c
A_bar = [0 1 0 0; 0 0 m*g/M 0; 0 0 0 1; 0 0 (g/(l*M))*(M + m) 0];
B_bar = [0; 1/M; 0; 1/(l*M)];
C_bar = [1 0 -l 0];
D_bar = [0];

% Construct A_hat matrix
size_A = size(A_bar);
size_B = size(B_bar);
zeros_A_hat = zeros(size_B(2), size_A(1) + size_B(2));
A_hat = [A_bar, B_bar; zeros_A_hat];

% Compute matrix exponential of A_hat matrix
matrix_exponential_A_hat = expm(A_hat * delta_t);

% Get F, G, H, M matrices
F = matrix_exponential_A_hat(1:4, 1:4);
G = matrix_exponential_A_hat(1:4, end);
H = C_bar;
M = D_bar;

% Stability
eigenvalues_F = eig(A_bar);

% Observability
O = [H; H*F; H*F*F; H*F*F*F];
observ = rank(O);


% Part d
% Construct L and y matrices as computed in part d
L = [H*(F - G*Kc); H*(F - G*Kc)^2; H*(F - G*Kc)^3; H*(F - G*Kc)^4];
y = [yNLhist(2); yNLhist(3); yNLhist(4); yNLhist(5)];

% Compute x_hat_0
x_hat_0 = inv(L' * L) * L' * y;

% Using x_hat_0, get first u_0 and y_0
u_0 = -Kc * x_hat_0;
y_0 = H * x_hat_0;

% Assign first approximations to predicted vectors
x_predicted(:,1) = x_hat_0;
u_predicted(1) = u_0;
y_predicted(1) = y_0;

% Go through time history and calculate predicted vectors using CT model
for i = 2:length(thist)
    x_predicted(:,i) = F * x_predicted(:,i-1) + G * u_predicted(i-1);
    u_predicted(i) = -Kc * x_predicted(:,i);
    y_predicted(i) = H * x_predicted(:,i);
end

% Plot figure
figure()
plot(thist, yNLhist, 'LineWidth',1.2)
hold on
plot(thist, y_predicted, '--', 'LineWidth',1.2)
legend("Actual History", "Predicted")
title("System Sensor Output")
xlabel("Time [seconds]")
ylabel("Pendulum bob's horizontal displacement [m]")



