clear; clc; close all;

% ASEN 5050 - Fall 2024
% Midterm 2 - Problem 2
% Jash Bhalavat

mu_jupiter = 1.268e8; % km3/s2
g0 = 9.81; % m/s2

% s/c orbiting jupiter
m_i = 2500; % kg
m_p = 500; % kg
I_sp = 300; % s
mu_X = 4000; % km3/s2

% Part a
v_in_i = [2.2789, 5.8841, 0]; % km/s
a_in_i = 1.8e6; % km
v_in_i_norm = norm(v_in_i);
eps_in_i = -mu_jupiter/(2*a_in_i);
r_in_i = mu_jupiter/(v_in_i_norm^2/2 - eps_in_i);
r_moon_x = r_in_i;

% Part b
v_moon_x = [0, sqrt(mu_jupiter/r_moon_x), 0];
v_inf_in = v_in_i - v_moon_x;

% Part c
v_out_norm = 7.6759;
phi_out = deg2rad(20.91);
v_r_out = v_out_norm * sin(phi_out);
v_theta_out = v_out_norm * cos(phi_out);
v_out = [v_r_out, v_theta_out, 0];

% Part d
v_inf_out = v_out - v_moon_x;

% Part h
v_inf_mag = 2.75;
eps_h = norm(v_inf_mag)^2/2;
a_h = -mu_X/(2*eps_h);
delta_v_eq = v_out - v_in_i;
delta = 2*asin((norm(delta_v_eq))/(2*v_inf_mag));
e_h = 1/(sin(delta/2));
r_p_h = a_h * (1 - e_h);

% Part i
delta_v_eq_norm = norm(delta_v_eq);

% Part j
delta_v_eq_norm_mps = 1000 * delta_v_eq_norm;
m_prop_needed = m_i * (1 - exp((-delta_v_eq_norm_mps)/(I_sp*g0)));

