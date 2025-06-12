close all; clear; clc;

% ASEN 6060 - HW 2, Problem 4
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
syms lambda_1 lambda_3 uxx

alpha_1 = (lambda_1^2 - uxx)/(2*lambda_1);
alpha_3 = (lambda_3^2 - uxx)/(2*lambda_3);

mat = [1 1 1 1;
        lambda_1 -lambda_1 lambda_3 -lambda_3;
        alpha_1 -alpha_1 alpha_3 -alpha_3;
        alpha_1*lambda_1 alpha_1*lambda_1 alpha_3*lambda_3 alpha_3*lambda_3];

mat_inv = inv(mat);

%% Part c

mu = mass_ratio_em;

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Calculate out of plane modes for all 5 eq points
for i = 1:5
    em_eq_pts_out_of_plane_modes(i,:) = out_of_plane_modes(mu, em_eq_pts(i,:));
    em_eq_pts_in_plane_modes(i,:) = in_plane_modes(mu, em_eq_pts(i,:));
end

% Only looking at L1 eq point planar oscillatory modes
l1_pos = em_eq_pts(1,:);
lambda_1 = em_eq_pts_in_plane_modes(1,1);
lambda_3 = em_eq_pts_in_plane_modes(1,3);
uxx_l1 = u_xx(mu, l1_pos);
uxy_l1 = u_xy(mu, l1_pos);
uyy_l1 = u_yy(mu, l1_pos);
alpha_1 = (lambda_1^2 - uxx_l1)/(2*lambda_1);
alpha_3 = (lambda_3^2 - uxx_l1)/(2*lambda_3);

xi_0 = -0.001;
% xi_0 = -0.000001;
eta_0 = 0.0;
xi_dot_0 = lambda_3 * eta_0 / alpha_3;
eta_dot_0 = alpha_3 * lambda_3 * xi_0;
init_var = [xi_0, xi_dot_0, eta_0, eta_dot_0];

% Time is one period
t = linspace(0, 2*pi/imag(lambda_3), 1000);

A3 = 1/(lambda_1^2 - lambda_3^2) * (xi_0*alpha_1*lambda_1 + (alpha_1*lambda_3*lambda_1*xi_dot_0)/uxx_l1 - (lambda_1^2*lambda_3*eta_0)/uxx_l1 - eta_dot_0);
A4 = 1/(lambda_1^2 - lambda_3^2) * (xi_0*alpha_1*lambda_1 - (alpha_1*lambda_3*lambda_1*xi_dot_0)/uxx_l1 + (lambda_1^2*lambda_3*eta_0)/uxx_l1 - eta_dot_0);

for i = 1:length(t)
    xi_t(i) = A3*exp(lambda_3*t(i)) + A4*exp(-lambda_3*t(i));
    eta_t(i) = A3*alpha_3*exp(lambda_3*t(i)) - A4*alpha_3*exp(-lambda_3*t(i));
end

figure()
plot(xi_t+l1_pos(1), eta_t, 'LineWidth',2)
hold on

x0 = [l1_pos(1) + xi_0; l1_pos(2) + eta_0; 0; xi_dot_0; eta_dot_0; 0];

% Set options for ode113
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Call ode113 function
[tout, xout] = ode113(@(t, state)CR3BP(state, mu), [0 t(end)], x0, options);

plot(xout(:,1), xout(:,2), 'LineWidth',2)
scatter(l1_pos(1), l1_pos(2), 'black', 'filled')
hold off
legend("Linearized Propagation", "CR3BP Propagation", "L1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
title("Periodic orbit around L1 in the linearized system and CR3BP")

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


%% Part e

U_star_XX = [uxx_l1, uxy_l1; uxy_l1, uyy_l1];
Omega = [0 2; -2 0];

A2D = [zeros(2), eye(2); U_star_XX, Omega];

[V, D] = eig(A2D);

evec_3 = V(:,3);
evec_4 = V(:,4);

basis = evec_3 + evec_4;


