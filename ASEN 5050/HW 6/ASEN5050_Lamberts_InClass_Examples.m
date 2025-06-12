clear; clc; close all;

% ASEN 5050, HW6
% Fall 2024
% Jash Bhalavat

% In class lambert's example

mu = 3.986004418e5;
R1 = [-654, 13605, 1997];
V1 = [-5.53, 0.849, 0.6830];
R2 = [7284, -19341, -3264];
V2 = [3.07, 2.63, 0.444];

TOF = 5*3600;
delta_a = 100000;
greater_than_180 = true;

delta_v_tot = lamberts_problem_solver(mu, R1, V1, R2, V2, TOF, delta_a, greater_than_180);

