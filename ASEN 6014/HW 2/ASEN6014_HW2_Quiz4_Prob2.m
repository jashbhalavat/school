clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 4, Problem 2
% Jash Bhalavat

%% Given


m = 1; % kg
g = 9.81; % m/s2
k = 3;

time = 20; % sec

y_0 = 1;
y_dot_0 = 0;
e1_t0 = y_0;
e2_t0 = y_dot_0;
e_t0 = [e1_t0; e2_t0];

[tout, state_out] = ode45(@(t, state)diff_eq(t, state, g, k, m), [0, time], e_t0);

function state_dot = diff_eq(t, state, g, k, m)
    % state - [e1, e2]
    omega = sqrt(k/m);
    P_12 = (cos(omega*t))^2 - (sin(omega*t))^2;
    P_21 = -P_12;
    P = [0, P_12; P_21, 0];

    dRde = [-g*cos(omega*t); -g/omega*sin(omega*t)];
    state_dot = P' * dRde;
end

% function state_dot = diff_eq(t, state, g, k, m)
%     omega = sqrt(k/m);
%     P_12 = 1/(cos(omega*t))^2 - 1/(sin(omega*t))^2;
%     P_21 = -P_12;
%     P = [0, P_12; P_21, 0];
% 
%     dRde = [-g*cos(omega*t); -g/omega * sin(omega*t)];
%     state_dot = P' * dRde;
% end


