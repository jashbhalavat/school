clear; clc; close all;

% ASEN 5050, HW 5, Problem 1
% Fall 2024, 10/17/2024
% Jash Bhalavat

mu_moon = 4902.799;
r_eq_moon = 1738;

% Lunar orbit s/c
a1 = 7045;
e1 = 0.23;

% Maneuver applied at -142 deg TA
theta_star_1 = deg2rad(-142);

% Part a
h1 = sqrt(mu_moon * (a1 * (1 - e1^2)));
v_r_1 = mu_moon/h1 * e1 * sin(theta_star_1);
v_theta_1 = mu_moon/h1 * (1 + e1 * cos(theta_star_1));
v1_rth = [v_r_1, v_theta_1, 0];
v1 = norm(v1_rth);

% Part b
delta_v = [0.3, -0.1, 0];

% Part c
v2_rth = delta_v + v1_rth;

% Part d
% r1 = r2
r1 = (h1^2/mu_moon) / (1 + e1*cos(theta_star_1));
r2 = r1;
v2 = norm(v2_rth);
eps2 = v2^2/2 - mu_moon/r2;
a2 = -mu_moon/(2*eps2);

% r = h^2/mu / (1 + ecos(theta_star))
% v_theta = mu/h *  (1 + ecos(theta_star))
% (1 + ecos(theta_star)) = v_theta * h/mu
% r = h^2/mu / v_theta*h/mu
% h = r * v_theta
v_theta_2 = v2_rth(2);
h2 = r2 * v_theta_2;
e2 = sqrt(1 + (2*h2^2*eps2)/mu_moon^2);

theta_star_2 = abs(acos(1/e2 * ((h2^2/mu_moon)/r2 - 1))) * sign(v2_rth(1));

% Part e
phi_fpa_1 = atan(v1_rth(1)/v1_rth(2));
phi_fpa_2 = atan(v2_rth(1)/v2_rth(2));
delta_phi_fpa = phi_fpa_2 - phi_fpa_1;
delta_v_mag = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(delta_phi_fpa));

% delta_v_mag = norm(delta_v) - this means the maneuver is coplanar
% i.e. change in theta_star = change in argument of periapsis
delta_theta_star = 2*pi - abs(theta_star_2) - abs(theta_star_1);

% Part f
m_tot = 1224;
m_prop = 156;
Isp = 212;
g = 9.81;
delta_v_norm_mps = norm(delta_v) * 1000;

% Using the rocket equation
m_prop_needed = m_tot * (1 - exp(-delta_v_norm_mps/(Isp * g)));



