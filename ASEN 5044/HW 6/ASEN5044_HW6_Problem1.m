clear; clc; close all;

% ASEN 5044, Fall 2024
% HW 6, Problem 1
% Jash Bhalavat

A = [0 1; -100 -10];
gamma = [0; 1];
W = 10;

n = length(A);

% Part b
delta_T = 0.2;

F = expm(A*delta_T);
zeros_z = zeros([n, n]);

Z = [-A, gamma * W * gamma'; zeros_z, A'] * delta_T;

e_z = expm(Z);

F_T = e_z(3:4, 3:4);
F_inv_Q = e_z(1:2, 3:4);

Q = F_T' * F_inv_Q;


