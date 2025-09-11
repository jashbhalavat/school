clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 5, Problem 1
% CWH Equations
% Jash Bhalavat

%% Given

a = 6800; % km

r_d_H = [1.299038,-1.0000,0.3213938]'; % km
v_d_H = [-0.000844437,-0.00292521, -0.000431250]'; % km/s

state_0 = [r_d_H; v_d_H];

time = 1300; % sec

[tout, state_out] = ode45(@(t,state)cwh_eqns(t,state,a), [0, time], state_0);

% # adjust the return matrix values as needed
% def result():
%     rho_H = [-0.606831046423625, -2.24354034342486, -0.346469112244555]  # km
%     rhoPrime_H = [-0.00154449565589346, 0.00136647986403504, -0.000405890168483601] # km/s
%     return rho_H, rhoPrime_H

function state_dot = cwh_eqns(t, state, a)
    % state - [r_d_H, v_d_H]
    x = state(1);
    y = state(2);
    z = state(3);
    xdot = state(4);
    ydot = state(5);
    zdot = state(6);

    mu = 398600;
    n = sqrt(mu/a^3);

    xdotdot = 2*n*ydot + 3*n^2*x;
    ydotdot = -2*n*xdot;
    zdotdot = -n^2*z;

    state_dot(1,1) = xdot;
    state_dot(2,1) = ydot;
    state_dot(3,1) = zdot;
    state_dot(4,1) = xdotdot;
    state_dot(5,1) = ydotdot;
    state_dot(6,1) = zdotdot;

end
