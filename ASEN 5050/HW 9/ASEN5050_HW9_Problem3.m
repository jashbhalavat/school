clear; clc; close all;

% ASEN 5050 - HW 9 - Problem 3
% Fall 2024
% Jash Bhalavat

mu_mars = 4.305e4; % km^3/s^2
P_rot_mars = 1.02595675*86400; % sec

omega_M = 2*pi/P_rot_mars;

lambda_omega = deg2rad(150);
P = lambda_omega/omega_M;

a = (mu_mars * (P/(2*pi))^2)^(1/3);
