clear; clc; close all;

% ASEN 5050, HW 5, Problem 2
% Fall 2024, 10/17/2024
% Jash Bhalavat

mu_mars = 4.305e4;
mu_moon = 4902.799;
mu_sun = 1.32712428e11;
mu_saturn = 3.794e7;

r_eq_mars = 3397.2;
r_eq_moon = 1738;

AU = 149597870.7;

a_earth = 1.0000010178 * AU;
a_saturn = 9.554909595 * AU;

% Given
% Mars orbit
% At t1 - before maneuver
r1 = 6500;
E1 = pi/2;
r_p_1 = 5915;

% At t2 - after maneuver
% s/c moving from periapsis to apoapsis
r_p_2 = 5712;
r_a_2 = 7888;

% Part a
% At t1, spacecraft is at b, r_b = a
a1 = r1;
e1 = 1 - r_p_1/a1;
p1 = a1 * (1 - e1^2);
h1 = sqrt(mu_mars * p1);

% Since E1 > 0, theta_star_1 > 0
theta_star_1 = abs(acos(1/e1 * ((h1^2/mu_mars)/r1 - 1)));

v_r_1 = mu_mars/h1 * e1 * sin(theta_star_1);
v_theta_1 = mu_mars/h1 * (1 + e1*cos(theta_star_1));

v1_rth = [v_r_1, v_theta_1, 0];


% Part b
r2 = r1;
e2 = (r_a_2 - r_p_2) / (r_a_2 + r_p_2);
a2 = r_a_2 / (1 + e2);
p2 = a2 * (1 - e2^2);
h2 = sqrt(mu_mars * p2);

% s/c from periapsis to apoapsis - theta_star > 0
theta_star_2 = abs(acos(1/e2 * ((h2^2/mu_mars)/r2 - 1)));

v_r_2 = mu_mars/h2 * e2 * sin(theta_star_2);
v_theta_2 = mu_mars/h2 * (1 + e2*cos(theta_star_2));

v2_rth = [v_r_2, v_theta_2, 0];

delta_v_rth = v2_rth - v1_rth;

delta_v_norm = norm(delta_v_rth);

v1 = norm(v1_rth);
v2 = norm(v2_rth);
phi_fpa_1 = atan(v1_rth(1)/v1_rth(2));
phi_fpa_2 = atan(v2_rth(1)/v2_rth(2));
delta_phi_fpa = phi_fpa_2 - phi_fpa_1;
delta_v_mag = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(delta_phi_fpa));


