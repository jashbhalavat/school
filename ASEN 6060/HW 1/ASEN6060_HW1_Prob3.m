clear; clc; close all;

% ASEN 6060 - HW 1, Problem 3
% Jash Bhalavat
% 01/28/2025

% Constants
G = 6.67408 * 10^-11; % m3/(kgs2)
G = G / (10^9); % km3/(kgs2)
% Earth
mu_earth = 398600.435507; % km3/s2
mass_earth = mu_earth / G; % kg
% Moon
mu_moon = 4902.800118; % km3/s2
a_moon = 384400; % km
mass_moon = mu_moon / G; % kg
% Earth-Moon system
mass_ratio_em = mass_moon / (mass_earth + mass_moon);
m_star_em = mass_earth + mass_moon;
mu = mass_ratio_em;

int_time = 25;
state0 = [0.12, 0, 0, 0, 3.45, 0];

% Set options for ODE45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) eventFn(t, y));

% Call ODE45 function
[tout, xout] = ode45(@(t, state)CR3BP(state, mu), [0 int_time], state0, options);

function state_dot = CR3BP(state, mu)
    % Circular Restricted 3 Body Problem non-dimensional EOMs
    x = state(1);
    y = state(2);
    z = state(3);
    xdot = state(4);
    ydot = state(5);
    zdot = state(6);

    r1 = sqrt((x + mu)^2 + (y)^2 + (z)^2);
    r2 = sqrt((x - 1 + mu)^2 + (y)^2 + (z)^2);

    state_dot(1, 1) = xdot;
    state_dot(2, 1) = ydot;
    state_dot(3, 1) = zdot;

    state_dot(4, 1) = 2*ydot + x - (1 - mu)*(x + mu)/(r1^3) - mu * (x - 1 + mu)/(r2^3);
    state_dot(5, 1) = -2*xdot + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3);
    state_dot(6, 1) = - (1 - mu)*z/(r1^3) - mu*z/(r2^3);
end

function [value,isterminal,direction] = eventFn(t,y)
    value = y(2); % Want y to be 0
    isterminal = 1; % Halt integration when value is 0
    direction = 1; % When zero is approached from +ve i.e. y_dot > 0
end

% Plot
figure()
plot3(xout(:,1), xout(:,2), xout(:,3))
hold on
scatter3(xout(1,1), xout(1,2), xout(1,3), 'filled')
scatter3(xout(end,1), xout(end,2), xout(end,3), 'filled', 'blue')
hold off
legend("Trajectory", "Initial State", "Final State")
title("Initial Condition - ", num2str(state0) )
xlabel("x")
ylabel("y")
zlabel("z")