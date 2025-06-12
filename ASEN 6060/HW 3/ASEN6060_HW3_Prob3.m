clear; clc; close all;

% ASEN 6060 - HW 3 Problem 3
% Spring 2025
% Jash Bhalavat

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
dx0 = [xi_0; eta_0; xi_dot_0; eta_dot_0];

% Set options for ode113
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

x0 = [l1_pos(1) + xi_0; l1_pos(2) + eta_0; 0; xi_dot_0; eta_dot_0; 0];

% Call ode113 function
[tout, xout] = ode113(@(t, state)CR3BP(state, mu), [0 t(end)], x0, options);


%% Part b

V0 = [x0; t(end)];
statef_V0 = xout(end,:);

V = V0;
Vi = V(:,1);
x0 = Vi(1:6);
phi0 = reshape(eye(6), [36,1]);
state0 = [x0; phi0];

[t_out, state_out] = ode113(@(t,state)CR3BP_full(state, mu), [0, V(7,1)], state0, options);
statef = state_out(end,:)';
state_dot = CR3BP_full(statef, mu);
statef_Vi = [statef(1:6); t_out(end)];

V_soln = gen_3d_periodic_orbit_single_shooting(Vi, mu, false);
[tout_corrected, xout_corrected] = ode113(@(t,state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

V_family = natural_param_continuation(V_soln, mu);


%% Plot all family members

load("prob3_tol_5e14.mat")

p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure()
scatter(l1_pos(1), l1_pos(2), 'filled', 'red')
hold on
scatter(p1_pos(1), p1_pos(2), 'filled', 'blue')
scatter(p2_pos(1), p2_pos(2), 'filled', 'black')

for i = 1:10:size(V_family, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family(7,i)], V_family(1:6,i), options);
    plot(xout(:,1), xout(:,2), 'LineWidth',2)
end
hold off
legend("L1", "Earth", "Moon")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
grid on
axis equal
title("L1 Lyapunov Orbit Family")

figure()
for i = 1:size(V_family,2)
    x = V_family(1,i);
    y = V_family(2,i);
    z = V_family(3,i);
    x_dot = V_family(4,i);
    y_dot = V_family(5,i);
    z_dot = V_family(6,i);
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
    C(i) = x^2 + y^2 + z^2 + 2*(1-mu)/r1 + mu/r2 - (x_dot^2 + y_dot^2 + z_dot^2);
end
plot(V_family(end,:)*t_star_em/86400, C, 'o')
hold on
scatter(V_family(end,1)*t_star_em/86400, C(1), 'filled', 'red')
scatter(V_family(end,end)*t_star_em/86400, C(end), 'filled', 'black')
legend("Jacobi Constant", "Start", "End")
xlabel("Period [days]")
ylabel("Jacobi Constant [-]")
grid on
title("Jacobi Constant and Period along the Family")
