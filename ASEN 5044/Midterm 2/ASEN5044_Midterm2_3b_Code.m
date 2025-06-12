clear; clc; close all;

% Given
delta_t = 0.5;
omega_a = 0.045;
H = [1 0 0 0; 0 0 1 0];
R_a_static = [75 7.5; 7.5 75];

% Helper variable for F_a matrix
odt_a = delta_t*omega_a;

% F_a matrix
F_a = [1 sin(odt_a)/omega_a 0 -(1-cos(odt_a))/omega_a;
        0 cos(odt_a) 0 -sin(odt_a);
        0 (1-cos(odt_a))/omega_a 1 sin(odt_a)/omega_a;
        0 sin(odt_a) 0 cos(odt_a)];

% Load data
data_part_b = load("midterm2_problem3b.mat");
y_a = data_part_b.yaHist;

% Get p, and T
[p, T] = size(y_a);

% Convert given measurements to column vector
y_vec = reshape(y_a, [T*p, 1]);

for i = 1:T
    % R_a is dynamically changing
    R_a_k((2*i-1):(2*i),:) = R_a_static + [12.5*sin(i/10), 25.5*sin(i/10); 25.5*sin(i/10), 12.5*cos(i/10)];
    
    % Assign temp variable
    R_a_k_temp = R_a_k((2*i-1):(2*i),:);

    % Assign H matrix
    H_mat((2*i-1):(2*i),:) = H*F_a^i;

    % Assign R matrix
    R_mat((2*i-1):(2*i), (2*i-1):(2*i)) = R_a_k_temp;
end

% Invert R matrix
R_mat_inv = inv(R_mat);

% State Estimate error covariance
P = inv(H_mat' * R_mat_inv * H_mat);

% Compute state initial estimate
x_hat_0 = P * H_mat' * R_mat_inv * y_vec;


