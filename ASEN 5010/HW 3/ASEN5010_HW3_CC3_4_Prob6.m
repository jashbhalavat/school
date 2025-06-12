clear; clc; close all;

% ASEN 5010 - HW 3, Concept Check 3, 4 Davenport's Q-Method, Problem 6
% Spring 2025
% Jash Bhalavat

v1_B = [0.8273; 0.5541; -0.0920];
v1_B_hat = normalize(v1_B);
v2_B = [-0.8285; 0.5522; -0.0955];
v2_B_hat = normalize(v2_B);

v1_N = [-0.1517; -0.9669; 0.2050];
v1_N_hat = normalize(v1_N);
v2_N = [-0.8393; 0.4494; -0.3044];
v2_N_hat = normalize(v2_N);

B = v1_B_hat * v1_N_hat' + v2_B_hat * v2_N_hat';
Z = [B(2,3) - B(3,2); B(3,1) - B(1,3); B(1,2) - B(2,1)];
S = B + B';
sigma = trace(B);

K = [sigma, Z';
    Z, S - sigma*eye(3)];

[V, D] = eig(K);

max_evec = V(:,4);

BN_est = ep2dcm(max_evec);

% dcm = [0.415936375258137, -0.854893534518477, 0.310087046449262,
%        -0.833756672976915, -0.494636774164691, -0.245324829380180,
%        0.363107066859884, -0.156497623887562, -0.918510616005040]

function hat = normalize(vec)
    vec_norm = norm(vec);
    hat = vec ./ vec_norm;
end

function out = ep2dcm(beta)
    beta0 = beta(1);
    beta1 = beta(2);
    beta2 = beta(3);
    beta3 = beta(4);
    out = [beta0^2 + beta1^2 - beta2^2 - beta3^2, 2*(beta1*beta2 + beta0*beta3), 2*(beta1*beta3 - beta0*beta2);
        2*(beta1*beta2 - beta0*beta3), beta0^2 - beta1^2 + beta2^2 - beta3^2, 2*(beta2*beta3 + beta0*beta1);
        2*(beta1*beta3 + beta0*beta2), 2*(beta2*beta3 - beta0*beta1), beta0^2 - beta1^2 - beta2^2 + beta3^2];
end