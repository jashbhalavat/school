clear; clc; close all;

% ASEN 6060 - HW 3 Problem 4
% Spring 2025
% Jash Bhalavat

%% Constants

G = 6.67408 * 10^-11; % m3/(kgs2)
G = G / (10^9); % km3/(kgs2)

% Earth
mu_earth = 398600.435507; % km3/s2
a_earth = 149598023; % km
e_earth = 0.016708617;
mass_earth = mu_earth / G; % kg

% Moon
mu_moon = 4902.800118; % km3/s2
a_moon = 384400; % km
e_moon = 0.05490;
mass_moon = mu_moon / G; % kg

% Earth-Moon system
mass_ratio_em = mass_moon / (mass_earth + mass_moon);
m_star_em = mass_earth + mass_moon;
l_star_em = a_moon;
t_star_em = sqrt(l_star_em^3/(G * m_star_em));

%% Correction

% Given
x0 = [0.82340, 0, -0.026755,0,0.13742,0]';
T = 2.7477;
V0 = [x0; T];

mu = mass_ratio_em;

% Set options for ode113
options = odeset('RelTol', 5e-14, 'AbsTol', 5e-14);

V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, true);

[tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Only looking at L1 eq point planar oscillatory modes
l1_pos = [em_eq_pts(1,:), 0];

% Set options for ode113
[tout, xout] = ode113(@(t, state)CR3BP(state, mu), [0 T], x0, options);

p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure()
hold on
scatter3(l1_pos(1), l1_pos(2), l1_pos(3), 'filled', 'red')
hold on
scatter3(x0(1), x0(2), x0(3), 'filled', 'blue')
plot3(xout_corrected(:,1), xout_corrected(:,2), xout_corrected(:,3), 'LineWidth',2)
% scatter3(p1_pos(1), p1_pos(2), p1_pos(3), 'filled', 'blue')
% scatter3(p2_pos(1), p2_pos(2), p2_pos(3), 'filled', 'black')
hold off
legend("L1", "Provided Initial Guess", "Corrected Trajectory")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)
axis equal
grid on
title("L1 Corrected Halo Orbit w/o the Moon")

figure()
plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2)
hold on
scatter3(l1_pos(1), l1_pos(2), l1_pos(3), 'filled', 'red')
scatter3(x0(1), x0(2), x0(3), 'filled', 'blue')
plot3(xout_corrected(:,1), xout_corrected(:,2), xout_corrected(:,3), 'LineWidth',2)
% scatter3(p1_pos(1), p1_pos(2), p1_pos(3), 'filled', 'blue')
scatter3(p2_pos(1), p2_pos(2), p2_pos(3), 'filled', 'black')
hold off
legend("CR3BP Propagation", "L1", "Initial Guess", "Corrected Trajectory")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)
axis equal
grid on
title("L1 Halo Orbit corrected trajectory")

