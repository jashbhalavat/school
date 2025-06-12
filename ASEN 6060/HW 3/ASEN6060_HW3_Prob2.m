clear; clc; close all;

% ASEN 6060 - HW 3 Problem 2
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

%% Part b

mu = mass_ratio_em;

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Only looking at L1 eq point planar oscillatory modes
l1_pos = [em_eq_pts(1,:), 0];

l1_in_plane_modes = in_plane_modes(mu, l1_pos);

oscillatory_eval = l1_in_plane_modes(3);
uxx_l1 = u_xx(mu, l1_pos);
uxy_l1 = u_xy(mu, l1_pos);
uyy_l1 = u_yy(mu, l1_pos);
U_star_XX = [uxx_l1, uxy_l1; uxy_l1, uyy_l1];

Omega = [0 2; -2 0];

A2D = [zeros(2), eye(2); U_star_XX, Omega];

[V, D] = eig(A2D);

oscillatory_evec = real(V(:,3));

oscillatory_pos_mag = norm([oscillatory_evec(1), oscillatory_evec(2)]);

pos_mag_req = 0.0001;

oscillatory_mag_factor = pos_mag_req / oscillatory_pos_mag;

oscillatory_ic = oscillatory_evec .* oscillatory_mag_factor;

% Time is one period
t = linspace(0, 2*pi/imag(oscillatory_eval), 1000);

xi_0 = oscillatory_ic(1);
xi_dot_0 = oscillatory_ic(3);
eta_0 = oscillatory_ic(2);
eta_dot_0 = oscillatory_ic(4);
x0 = [l1_pos(1) + xi_0; l1_pos(2) + eta_0; 0; xi_dot_0; eta_dot_0; 0];

% Set options for ode113
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[tout, xout] = ode113(@(t, state)CR3BP(state, mu), [0 t(end)], x0, options);

%% Part c
V0 = [x0; t(end)];

V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, true);

[tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

figure()
scatter(l1_pos(1), l1_pos(2), 'filled', 'black')
hold on
scatter(V_soln(1), V_soln(2), 'filled', 'blue')
plot(xout_corrected(:,1), xout_corrected(:,2), 'LineWidth',2)
hold off
axis equal
legend("L1", "Initial State", "Corrected Trajectory")
grid on
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
hold off
title("L1 Lyapunov Orbit - Corrected Trajectory")



figure()
plot(xout(:,1), xout(:,2), 'LineWidth',2)
hold on
scatter(l1_pos(1), l1_pos(2), 'filled', 'black')
scatter(xi_0 + l1_pos(1), eta_0 + l1_pos(2), 'filled', 'blue')

plot(xout_corrected(:,1), xout_corrected(:,2), 'LineWidth',2)
axis equal
legend("CR3BP Propagation", "L1", "Initial State", "Corrected Trajectory")
grid on
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
hold off
title("L1 Lyapunov Orbit - Initial and Corrected Trajectory")