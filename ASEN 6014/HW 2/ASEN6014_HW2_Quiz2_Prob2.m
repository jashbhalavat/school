clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 2, Problem 2
% Jash Bhalavat

%% Given

m = 1; % kg
g = 9.81; % m/s2
k = 3;

time = 10; % sec

y_0 = 1;
y_dot_0 = 0;
e1_t0 = y_0;
e2_t0 = y_dot_0;
e_t0 = [e1_t0; e2_t0];

[tout, state_out] = ode45(@(t, state)diff_eq(t, state, g, k, m), [0, time], e_t0);

% # adjust the return matrix values as needed
% def result():
%     L12 = [1.0]  
%     e1 = [4.13349747217459]  
%     e2 = [5.65845329933522]  
%     return L12, e1, e2




function state_dot = diff_eq(t, state, g, k, m)
    % state - [e1, e2]
    % L_12 = (cos(sqrt(k/m)*t))^2 + sqrt(k/m)*(sin(sqrt(k/m)*t))^2;
    L_12 = 1.0;
    L_21 = -L_12;
    L = [0, L_12; L_21, 0];

    % y_t = state(1)*cos(sqrt(k/m)*t) + state(2)*sin(sqrt(k/m)*t);
    % ydot_t = -sqrt(k/m)*state(1)*sin(sqrt(k/m)*t) + state(2)*cos(sqrt(k/m)*t);
    % dRde = [-g*y_t; -g*ydot_t];
    dRde = [-g*cos(sqrt(k/m)*t); -g/sqrt(k/m)*sin(sqrt(k/m)*t)];
    state_dot = inv(L) * dRde;
end