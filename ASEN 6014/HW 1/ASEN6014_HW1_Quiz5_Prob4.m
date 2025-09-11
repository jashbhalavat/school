clear; clc; close all;

% ASEN 6014 - HW 1, Quiz 5, Problem 4
% Fall 2025
% Jash Bhalavat

%% Given

r1_t0 = [-6685.20926,601.51244,3346.06634]; % km
v1_t0 = [-1.74294,-6.70242,-2.27739]; % km/s
state1_0 = [r1_t0';  v1_t0'];

r2_t0 = [-6685.21657,592.52839,3345.6716]; % km
v2_t0 = [-1.74283,-6.70475,-2.27334]; % km/s
state2_0 = [r2_t0';  v2_t0'];

time = 4848; % sec

state_0 = [state1_0; state2_0];
[tout, state_out] = ode45(@(t, state)dual_sc_diff_eq_J2(t, state), [0, time], state_0);

state_f = state_out(end,:);

function state_dot = dual_sc_diff_eq_J2(t, state)
    mu = 398600; % km3/s2
    J2 = 1082.63e-6;
    req = 6378.0; % km

    r1_N = state(1:3);
    v1_N = state(4:6);
    r2_N = state(7:9);
    v2_N = state(10:12);

    r1 = norm(r1_N);
    r2 = norm(r2_N);
    aj2_1 = -3/2 * J2 * mu/r1^2 * (req/r1)^2 * (1 - 5*(r1_N(3)/r1)^2) * r1_N/r1;
    aj2_2 = -3/2 * J2 * mu/r2^2 * (req/r2)^2 * (1 - 5*(r2_N(3)/r2)^2) * r2_N/r2;

    r1_dotdot = -mu/r1^3 * r1_N + aj2_1;
    r2_dotdot = -mu/r2^3 * r2_N + aj2_2;

    state_dot(1:3,1) = v1_N;
    state_dot(4:6,1) = r1_dotdot;
    state_dot(7:9,1) = v2_N;
    state_dot(10:12,1) = r2_dotdot;
end

% # adjust the return matrix values as needed
% def result():
%     r1 = [1750.36229893374, 6897.68437032018, 2363.64701098761]  # km
%     v1 = [-6.509533954556, 0.541555453028493, 3.23786390759598]  # km/sec
%     r2 = [1758.09754830167, 6889.25365625185, 2350.93087073889]  # km
%     v2 = [-6.5141781259016, 0.552148185898044, 3.2485081442567]  # km/sec
% 
%     return r1, v1, r2, v2
