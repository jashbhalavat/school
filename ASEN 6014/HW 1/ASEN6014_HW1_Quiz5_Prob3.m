clear; clc; close all;

% ASEN 6014 - HW 1, Quiz 5, Problem 3
% Fall 2025
% Jash Bhalavat

%% Given

r_t0 = [2466.69,5941.54,3282.71]; % km
v_t0 = [-6.80822,1.04998,3.61939]; % km/s
state_0 = [r_t0';  v_t0'];

time = 60 * 60;  % seconds

[tout, state_out] = ode45(@(t, state)diff_eq(t, state), [0, time], state_0);

r_f = state_out(end,1:3);
v_f = state_out(end,4:6);

function state_dot = diff_eq(t, state)
    mu = 398600; % km3/s2
    r_N = state(1:3);
    r_N_dot = state(4:6);

    r = norm(r_N);

    r_N_dotdot = -mu/r^3 * r_N;

    state_dot(1:3, 1) = r_N_dot;
    state_dot(4:6, 1) = r_N_dotdot;
end



% # adjust the return matrix values as needed
% def result():
%     rVec_N = [-3436.21920760318, -7152.73335600499, -3756.09636654493]  # km
%     vVec_N = [5.57477962735197, -0.921505143154174, -3.00853976950512]  # km/sec
% 
%     return rVec_N, vVec_N


