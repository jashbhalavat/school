clear; clc; close all;

% ASEN 5050 - HW 7
% Fall 2024
% Problem 1


mu_sun = 1.32712428e11;
mu_earth = 3.986004415e5;
mu_saturn = 3.794e7;
mu_jupiter = 1.268e8;

AU = 149597870.7;

a_earth = 1.0000010178 * AU;
a_saturn = 9.554909595 * AU;
a_jupiter = 5.202603191 * AU;

r_saturn = 60268;
r_jupiter = 71492;

% Consider a spacecraft in a large orbit around Saturn, described by a periapsis radius of 600,000
% km and an apoapsis radius of 1,800,000 km – and lying in the same orbit plane as Titan. Assume
% that Titan is in a circular orbit of radius 1,221,830 km, possesses a mass of 1.3455x1023 kg and is
% modeled as a sphere with a radius equal to 2,575 km. Also assume that the plane of the
% spacecraft’s orbit relative to Saturn does not change through the flyby.

r_p = 600000;
r_a = 1800000;

titan_orbit_radius = 1221830;
titan_mass = 1.3455e23;
titan_radius = 2575;

G = 6.67430e-20;
mu_titan = G * titan_mass;

%% Part a
e = (r_a - r_p) / (r_a + r_p);
a = 1/2 * (r_a + r_p);
p = a * (1 - e^2);
h = sqrt(p * mu_saturn);
b = a * sqrt(1-e^2);

r_intersection_w_titan = titan_orbit_radius;
theta_star = -acos(1/e * (p/r_intersection_w_titan - 1));

v_r = mu_saturn/h * e * sin(theta_star);
v_theta = mu_saturn/h * (1 + e * cos(theta_star));
v_sc_intersection = [v_r, v_theta, 0];

v_titan_theta = sqrt((mu_saturn) / titan_orbit_radius);
v_titan = [0 v_titan_theta 0];
v_titan_mag = norm(v_titan);

v_inf_in = v_sc_intersection - v_titan;

%% Part b
r_p_h = 3000;

v_inf_in_mag = norm(v_inf_in);

eps_h = v_inf_in_mag^2/2;
a_h = -mu_titan/(2*eps_h);

e_h = 1 - r_p_h/a_h;

delta = 2*asin(1/e_h);

%% Part e
beta_1 = pi - acos((dot(v_inf_in, v_titan))/(v_inf_in_mag*v_titan_mag));
beta_2 = beta_1 - delta;

v_inf_out_mag = v_inf_in_mag;

v_out_squared = v_inf_out_mag^2 + v_titan_mag^2 - 2*v_inf_out_mag*v_titan_mag*cos(beta_2);
v_out_norm = sqrt(v_out_squared);

phi_fpa = -acos((v_inf_out_mag^2 - v_out_norm^2 - v_titan_mag^2)/(-2*v_out_norm*v_titan_mag));

v_out_r = v_out_norm * sin(phi_fpa);
v_out_theta = v_out_norm * cos(phi_fpa);
v_out = [v_out_r, v_out_theta, 0];

%% Part f

h_out = r_intersection_w_titan * v_out_theta;
eps_out = v_out_norm^2/2 - mu_saturn/r_intersection_w_titan;
a_out = -mu_saturn/(2*eps_out);
e_out = sqrt((eps_out*2*h_out^2)/(mu_saturn^2) + 1);
theta_star_out = -acos(1/e_out * ((h_out^2/mu_saturn)/r_intersection_w_titan - 1));

% test_theta_star_out = asin((v_out_r*h_out)/(mu_saturn*e_out))
% test_theta_star_out_2 = acos((v_out_theta * h_out/mu_saturn - 1)/e_out)

%% Part g
delta_v_eq = v_out - v_sc_intersection;
norm_delta_v_eq = norm(delta_v_eq);

test_norm_delta_v_eq = 2 * v_inf_out_mag * sin(delta/2);










