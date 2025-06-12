clear; clc; close all;

% ASEN 5044, HW 7
% Fall 2024
% Jash Bhalavat

syms delta_t
syms omega

A = [0 1 0 0; 0 0 0 -omega; 0 0 0 1; 0 omega 0 0];

e_at = expm(A * delta_t);

delta_t = 0.5;
omega = 0.045;

% Omega times delta_t
odt = omega * delta_t;

e_at_given = [1 sin(odt)/omega 0 -(1 - cos(odt))/omega; 0 cos(odt) 0 -sin(odt); 0 (1-cos(odt))/omega 1 sin(odt)/omega; 0 sin(odt) 0 cos(odt)];

u_0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P_a_0 = diag([10, 2, 10, 2]);

P_1_2sig(1) = 2 * sqrt(P_a_0(1, 1));
P_2_2sig(1) = 2 * sqrt(P_a_0(2, 2));
P_3_2sig(1) = 2 * sqrt(P_a_0(3, 3));
P_4_2sig(1) = 2 * sqrt(P_a_0(4, 4));

T = 300;

F = e_at_given;

mu = [u_0];
P_k_minus_1 = [P_a_0];

for i = 1:T
    mu(:, i+1) = F * mu(:, i);
    P_k_minus_1 = F * P_k_minus_1 * F';
    P_1_2sig(i+1) = 2 * sqrt(P_k_minus_1(1, 1));
    P_2_2sig(i+1) = 2 * sqrt(P_k_minus_1(2, 2));
    P_3_2sig(i+1) = 2 * sqrt(P_k_minus_1(3, 3));
    P_4_2sig(i+1) = 2 * sqrt(P_k_minus_1(4, 4));
end

function out = get_diag(mat)
    len = length(mat);
    for i = 1:len
        out(i, 1) = mat(i, i);
    end
end


k = 0:T;

figure()
plot(k, mu(1, :))
hold on
plot(k, P_1_2sig)
plot(k, -P_1_2sig)
hold off


figure()
plot(k, mu(2, :))
hold on
plot(k, P_2_2sig)
plot(k, -P_2_2sig)
hold off

figure()
plot(k, P_1_2sig)

figure()
plot(k, P_2_2sig)
