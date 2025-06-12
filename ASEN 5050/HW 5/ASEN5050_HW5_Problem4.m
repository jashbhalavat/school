clear; clc; close all;

% ASEN 5050, HW 5, Problem 4
% Fall 2024, 10/17/2024
% Jash Bhalavat

mu_mars = 42828.3143;

a = 6600;
e = 0.46;

% Given
E1 = deg2rad(333.1636 - 360);
E2 = deg2rad(26.5500);

P = 2 * pi * sqrt(a^3/mu_mars);
n = (2*pi) / P;

t1_to_tp = abs(1/n * (E1 - e*sin(E1)));
tp_to_t2 = 1/n * (E2 - e*sin(E2));

total_t = t1_to_tp + tp_to_t2;


