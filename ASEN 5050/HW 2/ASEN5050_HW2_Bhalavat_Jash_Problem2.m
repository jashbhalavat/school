clear; clc; close all;

mu_moon = 4.902799e3;
r_eq_moon = 1738;

% Problem 2

n = [0.6428 -0.7660 0];
h = [-0.3237 -0.2717 0.9063];
h_norm = norm(h);
e = [0.0475 0.3755 0.1295];
e_norm = norm(e);


% Part a
e_norm = norm(e);

i = acos(h(3));

raan = abs(acos(n(1)));
raan = sign(n(2)) * raan;

aop = abs(acos(dot(n, e)/norm(e)));
aop = sign(e(3)) * aop;

% Part b
r = 4070.6;
theta_star = -aop;

h_an = sqrt(r * mu_moon * (1 + e_norm*cos(theta_star)));

a = (h_an^2/mu_moon)/(1 - e_norm^2);

