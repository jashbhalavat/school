clear; clc; close all;

% ASEN 6060 - HW 2, Problem 5
% Jash Bhalavat
% Spring 2025

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

%% Part a

mu = mass_ratio_em;

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Only looking at L4 eq point planar oscillatory modes
l4_pos = em_eq_pts(4,:);

l4_in_plane_modes = in_plane_modes(mu, l4_pos);

lambda_1 = l4_in_plane_modes(1);
lambda_3 = l4_in_plane_modes(3);
uxx_l4 = u_xx(mu, l4_pos);
uxy_l4 = u_xy(mu, l4_pos);
uyy_l4 = u_yy(mu, l4_pos);
U_star_XX = [uxx_l4, uxy_l4; uxy_l4, uyy_l4];

Omega = [0 2; -2 0];

A2D = [zeros(2), eye(2); U_star_XX, Omega];

[V, D] = eig(A2D);

alpha_1 = (lambda_1^2 - uxx_l4)/(2*lambda_1);
alpha_3 = (lambda_3^2 - uxx_l4)/(2*lambda_3);

% sp_evec = V(:,1) + V(:,2);
sp_evec = real(V(:,1));
lp_evec = V(:,3) + V(:,4);

% sp_evec_norm = sp_evec ./ norm(sp_evec);
% lp_evec_norm = lp_evec ./ norm(lp_evec);

% sp_pos_mag = norm([sp_evec_norm(1), sp_evec_norm(3)]);
% lp_pos_mag = norm([lp_evec_norm(1), lp_evec_norm(3)]);

sp_pos_mag = norm([sp_evec(1), sp_evec(2)]);
lp_pos_mag = norm([lp_evec(1), lp_evec(3)]);

pos_mag_req = 0.02;

sp_mag_factor = pos_mag_req / sp_pos_mag;
lp_mag_factor = pos_mag_req / lp_pos_mag;

sp_ic = sp_evec .* sp_mag_factor;
lp_ic = lp_evec .* lp_mag_factor;

% Short Period Linear Prop
xi_0 = sp_ic(1);
xi_dot_0 = sp_ic(3);
eta_0 = sp_ic(2);
eta_dot_0 = sp_ic(4);
% xi_dot_0 = lambda_3 * eta_0 / alpha_3;
% eta_dot_0 = alpha_3 * lambda_3 * xi_0;

% Time is one period
t = linspace(0, 2*pi/imag(lambda_3), 1000);

A3 = 1/(lambda_1^2 - lambda_3^2) * (xi_0*alpha_1*lambda_1 + (alpha_1*lambda_3*lambda_1*xi_dot_0)/uxx_l4 - (lambda_1^2*lambda_3*eta_0)/uxx_l4 - eta_dot_0);
A4 = 1/(lambda_1^2 - lambda_3^2) * (xi_0*alpha_1*lambda_1 - (alpha_1*lambda_3*lambda_1*xi_dot_0)/uxx_l4 + (lambda_1^2*lambda_3*eta_0)/uxx_l4 - eta_dot_0);

for i = 1:length(t)
    xi_t(i) = A3*exp(lambda_3*t(i)) + A4*exp(-lambda_3*t(i));
    eta_t(i) = A3*alpha_3*exp(lambda_3*t(i)) - A4*alpha_3*exp(-lambda_3*t(i));
end

figure()
plot(xi_t+l4_pos(1), eta_t+l4_pos(2))
hold on

x0 = [l4_pos(1) + xi_0; l4_pos(2) + eta_0; 0; xi_dot_0; eta_dot_0; 0];

% Set options for ODE45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Call ODE45 function
[tout, xout] = ode45(@(t, state)CR3BP(state, mu), [0 t(end)], x0, options);

plot(xout(:,1), xout(:,2))
plot(l4_pos(1), l4_pos(2))
hold off
legend("Linearized Propagation", "CR3BP Propagation", "L4")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
title("Short Periodic orbit around L4 in the linearized system and CR3BP")

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

% Long Period Linear Prop
xi_0 = lp_ic(1);
% xi_dot_0 = lp_ic(2);
eta_0 = lp_ic(3);
% eta_dot_0 = lp_ic(4);
xi_dot_0 = lambda_3 * eta_0 / alpha_3;
eta_dot_0 = alpha_3 * lambda_3 * xi_0;

% Time is one period
t = linspace(0, 2*pi/imag(lambda_3), 1000);

A3 = 1/(lambda_1^2 - lambda_3^2) * (xi_0*alpha_1*lambda_1 + (alpha_1*lambda_3*lambda_1*xi_dot_0)/uxx_l4 - (lambda_1^2*lambda_3*eta_0)/uxx_l4 - eta_dot_0);
A4 = 1/(lambda_1^2 - lambda_3^2) * (xi_0*alpha_1*lambda_1 - (alpha_1*lambda_3*lambda_1*xi_dot_0)/uxx_l4 + (lambda_1^2*lambda_3*eta_0)/uxx_l4 - eta_dot_0);

for i = 1:length(t)
    xi_t(i) = A3*exp(lambda_3*t(i)) + A4*exp(-lambda_3*t(i));
    eta_t(i) = A3*alpha_3*exp(lambda_3*t(i)) - A4*alpha_3*exp(-lambda_3*t(i));
end

figure()
plot(xi_t+l4_pos(1), eta_t+l4_pos(2))
hold on

x0 = [l4_pos(1) + xi_0; l4_pos(2) + eta_0; 0; xi_dot_0; eta_dot_0; 0];

% Set options for ODE45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Call ODE45 function
[tout, xout] = ode45(@(t, state)CR3BP(state, mu), [0 t(end)], x0, options);

plot(xout(:,1), xout(:,2))
plot(l4_pos(1), l4_pos(2))
hold off
legend("Linearized Propagation", "CR3BP Propagation", "L4")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
title("Long Periodic orbit around L4 in the linearized system and CR3BP")




