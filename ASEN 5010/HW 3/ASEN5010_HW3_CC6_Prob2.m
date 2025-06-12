clear; clc; close all;

% ASEN 5010 - HW 3, Concept Check 6 OLAE, Problem 2
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

s1 = v1_B_hat + v1_N_hat;
s2 = v2_B_hat + v2_N_hat;

d1 = v1_B_hat - v1_N_hat;
d2 = v2_B_hat - v2_N_hat;

S = [skew_symmetric(s1); skew_symmetric(s2)];
d = [d1; d2];
W = [eye(3), zeros(3);
    zeros(3), eye(3)];

q_bar = inv(S' * W * S) * S' * W * d;
beta_bar = (1/(sqrt(1 + q_bar' * q_bar))) * [1; q_bar];


BN_est = ep2dcm(beta_bar);

% dcm = [0.416218655892902, -0.854762273457525, 0.310070131358101,
%        -0.833607926580328, -0.494900912533443, -0.245297597860255,
%        0.363125123379190, -0.156379482820763, -0.918523599110855]

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

function out = skew_symmetric(vec)
    out = [0 -vec(3) vec(2);
            vec(3) 0 -vec(1);
            -vec(2) vec(1) 0];
end