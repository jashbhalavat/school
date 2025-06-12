clear; clc; close all;

% ASEN 5050 - HW 6
% Fall 2024, 10/31/2024
% Jash Bhalavat
% Problem 1

mu_sun = 1.32712428e11;
AU = 149597870.7;

% Earth to Venus
jd1 = 2457754.5;
jd2 = 2457871.5;

TOF = (jd2 - jd1) * 86400;

R1 = [-2.686982e7, 1.326980e8, 5.752566e7];
V1 = [-29.781722, -5.101774, -2.210394];

R2 = [-5.642901e7, -8.571048e7, -3.499466e7];
V2 = [29.658341, -16.091100, -9.116674];

delta_a = 0.1;
greater_than_180 = false;

delta_v_out = lamberts_problem_solver_ellipitical(mu_sun, R1, V1, R2, V2, TOF, delta_a, greater_than_180);
