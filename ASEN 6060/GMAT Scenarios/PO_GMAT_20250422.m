clear; clc; close all;

% ASEN 6060 - 04/22/2025 Periodic Orbit GMAT Scenario
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
mu = mass_ratio_em;

%%

x_bar = [8.35071e-1, 0.0, 1.40622e-1 , 0, 2.51487e-1 , 0];
T = 2.76312;

x0 = x_bar;
x0(1) = x0(1) - (1-mu);

R_moon_wrt_earth_GCRF = [1.890395e+04, 3.292575e+05, 1.785784e+05]; % km
l_star = norm(R_moon_wrt_earth_GCRF);
t_star = sqrt((l_star^3)/(G*m_star_em));
tau = 2*pi * sqrt(l_star^3/mu);

x0_dim(1:3) = x0(1:3) * l_star;
x0_dim(4:6) = x0(4:6) * l_star/t_star;
t_dim = tau * t_star;

% Set options for ODE45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

x0(1) = x0(1) + (1-mu);

[tout, xout] = ode45(@(t, state)CR3BP(state, mu), [0 T/2], x0, options);

x_after_half_P = xout(end,:);
x_after_half_P(1) = x_after_half_P(1) - (1-mu);

R_moon_wrt_earth_GCRF_at_half_P = [-3.33590e+05; 1.82242e+05; 9.83307e+04]; % km
l_star_half_P = norm(R_moon_wrt_earth_GCRF_at_half_P);
t_star_halp_P = sqrt((l_star_half_P^3)/(G*m_star_em));

x_after_half_P_dim(1:3) = x_after_half_P(1:3) * l_star_half_P;
x_after_half_P_dim(4:6) = x_after_half_P(4:6) * l_star_half_P/t_star_halp_P;

R_moon_wrt_earth_GCRF_at_P = [ -3.69926e+05; -1.42961e+05; -7.85727e+04];
l_star_P = norm(R_moon_wrt_earth_GCRF_at_P);
t_star_P = sqrt((l_star_P^3)/(G*m_star_em));

x0 = x_bar;
x0(1) = x0(1) - (1-mu);

x_after_P= x_bar;
x_after_P(1) = x_after_P(1) - (1-mu);

x_after_P_dim(1:3) = x_after_P(1:3) * l_star_P;
x_after_P_dim(4:6) = x_after_P(4:6) * l_star_P/t_star_P;

