clear; clc; close all;

mu_sun = 1.32712428e11;

r = 1.6599e8;
a = 2.3132e8;
v_r = -11.6485;

eps = -mu_sun/(2*a);
v = sqrt(mu_sun*(2/r - 1/a));
v_theta = sqrt(v^2 - v_r^2);

phi_fpa = atan(v_r/v_theta);
h = r*v*cos(phi_fpa);
p = h^2/mu_sun;
e = sqrt(1-p/a);
theta_star = -acos(((h^2/mu_sun)/r - 1)/e);
P = 2 * pi * sqrt(a^3/mu_sun);

b = a * sqrt(1 - e^2)

r_max = (h^2/mu_sun)/(1+e*cos(-pi));
r_min = (h^2/mu_sun)/(1+e*cos(0));
%%
clear; clc; close all;

mu_sun = 1.32712428e11;

R = [1.0751e8 -1.2647e8 1.3644e5];
V = [1.5180e1 2.8193e1 1.0504e-2];

% Part d
H = cross(R, V);
h = norm(H);

E = cross(V, H)/mu_sun - R/norm(R);
e = norm(E);

EPS = 1/2*norm(V)^2 - mu_sun/norm(R)