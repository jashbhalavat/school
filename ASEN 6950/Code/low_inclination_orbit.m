clear; clc; close all;

%% Constants and IC

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

mu = mass_ratio_em;

halo_orbits = load("V_family_halo_L1.mat").V_family;
em_eq_pts = load("eq_points.mat").em_eq_pts;

l2_pos = [em_eq_pts(2,:), 0];
l1_pos = [em_eq_pts(1,:), 0];
p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure(1)
scatter3(l1_pos(1), l1_pos(2), 0, 'filled', 'red')
hold on
scatter3(l2_pos(1), l2_pos(2), 0, 'filled', 'red')
scatter3(p1_pos(1), p1_pos(2), 0, 'filled', 'blue')
scatter3(p2_pos(1), p2_pos(2), 0, 'filled', 'black')

% Set options for ode113
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

n = 1;

[tout_halo, xout_halo] = ode113(@(t,state)CR3BP(state, mu), [0, halo_orbits(7,n)], halo_orbits(1:6,n), options);
plot3(xout_halo(:,1), xout_halo(:,2), xout_halo(:,3), 'LineWidth',2)

lyapunov_orbits = load("V_family_L2_Lyapunov_orbits.mat").V_family;
[tout_lyapunov, xout_lyapunov] = ode113(@(t,state)CR3BP(state, mu), [0, lyapunov_orbits(7,end-3)], lyapunov_orbits(1:6,end-3), options);
plot3(xout_lyapunov(:,1), xout_lyapunov(:,2), xout_lyapunov(:,3), 'LineWidth',2)

legend("L1", "L2", "Earth", "Moon", "Initial Halo Orbit")
title("Initial Halo Orbit (Low Inclination)")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)


