clear; clc; close all;

% ASEN 5050 - Fall 2024
% Midterm 2 - Problem 3
% Jash Bhalavat

mu_jupiter = 1.268e8; % km3/s2
mu_sun = 1.32712428e11; % km3/s2
AU = 149597870.7; % km

% Interplanetary mission, jupiter flyby
a_jupiter = 5.202603191 * AU;

v_in = [3.2476, 13.1175, 0];
v_out = [3.3426, 13.7356, 0];

v_theta_jupiter = sqrt(mu_sun/a_jupiter);
v_jupiter = [0, v_theta_jupiter, 0];

v_inf_in = v_in - v_jupiter;
v_inf_out = v_out - v_jupiter;

v_inf_in_norm = norm(v_inf_in);
v_inf_out_norm = norm(v_inf_out);

% Using v_inf_in
eps_h = v_inf_in_norm^2/2;
a_h = -mu_jupiter/(2*eps_h);
delta_v_eq = v_out - v_in;
delta = 2*asin((norm(delta_v_eq))/(2*v_inf_in_norm));
e_h = 1/(sin(delta/2));
r_p_h_in = a_h * (1 - e_h);

% Using v_inf_out
eps_h = v_inf_out_norm^2/2;
a_h = -mu_jupiter/(2*eps_h);
delta_v_eq = v_out - v_in;
delta = 2*asin((norm(delta_v_eq))/(2*v_inf_out_norm));
e_h = 1/(sin(delta/2));
r_p_h_out = a_h * (1 - e_h);

