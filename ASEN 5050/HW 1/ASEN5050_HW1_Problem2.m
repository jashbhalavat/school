clear; clc; close all;

v_inf = 10.7527;
theta_inf = deg2rad(139.3724);

mu_jupiter = 1.268e8;

a = - mu_jupiter / v_inf^2;
e = -1 / cos(theta_inf);
delta = 2 * asin(1/e);

r_p = abs(a) * (e^2 - 1) / (1 + e);
v_p = sqrt((v_inf^2/2 + mu_jupiter/r_p)*2);

v_test = sqrt((mu_jupiter/(2*abs(a)) + mu_jupiter/r_p)*2);
