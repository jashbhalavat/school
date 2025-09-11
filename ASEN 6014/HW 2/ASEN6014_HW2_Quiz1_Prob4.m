clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 1, Problem 4
% Jash Bhalavat

%% Given

t = 10;

e1_t0 = 1;
e2_t0 = 0;
k1 = 3;
k3 = 1;

e_t0 = [e1_t0; e2_t0];
k = [k1; k3];

[tout, state_out] = ode45(@(t, state)diff_eq(t, state, k), [0, t], e_t0);

function e_dot = diff_eq(t, state, k)
    % state - e1, e2
    % k - k1, k3
    x_t = state(1) * cos(sqrt(k(1))*t) + state(2)/sqrt(k(1)) * sin(sqrt(k(1))*t);
    ad = -k(2) * x_t^3;
    e1_dot = -(1/sqrt(k(1)) * sin(sqrt(k(1))*t)) * ad;
    e2_dot = cos(sqrt(k(1))*t) * ad;

    e_dot = [e1_dot; e2_dot];
end