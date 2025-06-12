clear; clc; close all;

% ASEN 5050 - HW 8, Problem 1
% Fall 2024
% Jash Bhalavat

AU = 149597870.7;
mu_sun = 1.32712428e11;
G = 6.67e-20;

r_p = 7500;
r_a = 8500;
i = deg2rad(105);
P_sc = 110 * 60;
r_planet = 6500;
a_planet = 2.25 * AU;

mars_orbit_period = 2*pi*sqrt(a_planet^3/mu_sun);
raan_dot = (360*pi/180)/mars_orbit_period;

a_sc = 1/2 * (r_p + r_a);
e_sc = (r_a - r_p)/(r_a + r_p);
mu_planet = (a_sc^3)/(P_sc/(2*pi))^2;
m_planet = mu_planet/G;

J_2 = (raan_dot * -2/3 * (1 - e_sc^2)^2 * a_sc^(7/2)) / (sqrt(mu_planet) * r_planet^2 * cos(i));


