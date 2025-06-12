clear; clc; close all;

% Given
mu_saturn = 3.794e7;
mars_r_eq = 3397.2;
saturn_r_eq = 60268;

% Problem 1
r1 = [-720000, 670000, 310000];
v1 = [2.160, -3.360, 0.620];

coe = rv2coe(r1, v1, mu_saturn);
a = coe(1);
e = coe(2);
i = coe(3);
raan = coe(4);
aop = coe(5);
theta_star = coe(6);

h = cross(r1, v1);
h_norm = norm(h);

% Point 2 is impact on surface
theta_star_2 = abs(acos(((h_norm^2/mu_saturn)/saturn_r_eq - 1)/e)) * -1;

p = a * (1 - e^2);
r2_pqw = [p*cos(theta_star_2)/(1 + e*cos(theta_star_2)), p*sin(theta_star_2)/(1 + e*cos(theta_star_2)), 0];
v2_pqw = [-sqrt(mu_saturn/p) * sin(theta_star_2), sqrt(mu_saturn/p)*(e + cos(theta_star_2)), 0];

r2_from_pqw = R3(-raan) * R1(-i) * R3(-aop) * r2_pqw';
v2_from_pqw = R3(-raan) * R1(-i) * R3(-aop) * v2_pqw';

coe(6) = theta_star_2;

rv = coe2rv(coe, mu_saturn);

r2_from_rth = rv(1:3);
v2_from_rth = rv(4:6);
