clear; clc; close all;

% ASEN 6060 - HW 1, Problem 2
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

init_state_case = 4;

switch init_state_case
    case 1
        % Integration time [unitless]
        int_time = 2;
        % Initial State [unitless]
        state0 = [0.98, 0, 0, 0, 1.2, 0];
    case 2
        int_time = 8;
        state0 = [0.98, 0, 0, 0, 1.7, 0];
    case 3
        int_time = 25;
        state0 = [0.12, 0, 0, 0, 3.45, 0];
    case 4
        int_time = 250;
        state0 = [0.12, 0, 0, 0, 3.48, 0];
end

% Set options for ODE45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

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

% Plot
figure()
plot3(xout(:,1), xout(:,2), xout(:,3))
hold on
scatter3(xout(1,1), xout(1,2), xout(1,3), 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 .75 .75])
scatter3(xout(end,1), xout(end,2), xout(end,3), 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.95 .95 .15])
scatter3(-mu, 0, 0, 'filled', 'blue')
scatter3(1-mu, 0, 0, 'filled', 'red')
hold off
legend("Trajectory", "Initial State", "Final State", "Earth", "Moon")
title("Initial State #" + init_state_case + " Trajectory")
xlabel("x")
ylabel("y")
zlabel("z")

%% Part d
int_time = 25;
state0 = [0.12, 0, 0, 0, 3.45, 0];
l_star = a_moon;
t_star = (l_star^3/(G*m_star_em))^(1/2);
initial_pos_dim = state0(1:3) * l_star;
initial_vel_dim = state0(4:6) * l_star/t_star;
tau = int_time * t_star;

% Planet pos
pos_p1 = mu * l_star;
sc_at_t0_from_earth = initial_pos_dim + [pos_p1, 0, 0];

t_terms_of_period = int_time / (2*pi);
