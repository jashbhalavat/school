clear; clc; close all;

% ASEN 5050 - Fall 2024
% Midterm 2 - Problem 1
% Jash Bhalavat

% Orbit around the moon
r_p_i = 1850; % km
r_a_i = 5400; % km
% Apoapsis to Periapsis
r_sc_1 = 3500; % km
r_sc_2 = 3500; % km
r_p_f = 3000; % km
r_a_f = 5400; % km

mu_moon = 4.902799e+03; % km3/s2

% Part a
% Initial orbital params
a_i = 1/2 * (r_p_i + r_a_i);
eps_i = -mu_moon/(2*a_i);
e_i = (r_a_i - r_p_i)/(r_a_i + r_p_i);
p_i = a_i * (1 - e_i^2);
h_i = sqrt(p_i * mu_moon);
theta_star_1 = -acos((p_i/r_sc_1 - 1)/e_i);

v_r_1 = mu_moon/h_i * e_i * sin(theta_star_1);
v_theta_1 = mu_moon/h_i * (1 + e_i*cos(theta_star_1));
v_1 = [v_r_1, v_theta_1, 0];


% Final orbital params
a_f = 1/2 * (r_p_f + r_a_f);
eps_f = -mu_moon/(2*a_f);
e_f = (r_a_f - r_p_f)/(r_a_f + r_p_f);
p_f = a_f * (1 - e_f^2);
h_f = sqrt(p_f * mu_moon);
theta_star_2 = -acos((p_f/r_sc_2 - 1)/e_f);

v_r_2 = mu_moon/h_f * e_f * sin(theta_star_2);
v_theta_2 = mu_moon/h_f * (1 + e_f*cos(theta_star_2));
v_2 = [v_r_2, v_theta_2, 0];

delta_v = v_2 - v_1;
delta_v_norm = norm(delta_v);

% Part b
change_in_aop = theta_star_2 - theta_star_1;





