clear; clc; close all;

mu_mars = 4.305e4;
r_eq_mars = 3397.2;

% Problem 1

% Part a
r = [3.62067e3 -3.19925e2 -4.20645e2];
r_norm = norm(r);
r_hat = r/r_norm;

v = [-4.28843e-1 -3.00176e-2 -3.39801];
v_norm = norm(v);
v_hat = v/v_norm;

h = cross(r, v);
h_norm = norm(h);
h_hat = h/h_norm;

K = [0 0 1];
n = cross(K, h);
n_norm = norm(n);
n_hat = n/n_norm;

e = 1/mu_mars * ((v_norm^2 - mu_mars/r_norm)*r - dot(r,v)*v);
e_norm = norm(e);
e_hat = e/e_norm;

eps = v_norm^2/2 - mu_mars/r_norm;

if e == 1.0
    p = h_norm^2/mu_mars;
else
    a = -mu_mars/(2*eps);
    p = a * (1 - e_norm^2);
end

i = acos(h_hat(3));
raan = sign(n_hat(2)) * abs(acos(n_hat(1)));
aop = sign(e(3)) * abs(acos(dot(n, e)/(n_norm * e_norm)));
theta_star = sign(dot(r, v)) * abs(acos(dot(e, r)/(e_norm * r_norm)));

% Special cases
aop_true = acos(e_hat(1));
if e(2) < 0
    aop_true = 2*pi - aop_true;
end


% Part b1
theta = theta_star + aop;

C_XYZ_RTH = R3(-raan) * R1(-i) * R3(-theta);
r_rth = C_XYZ_RTH' * r';
v_rth = C_XYZ_RTH' * v';

% Part b2
v_r = mu_mars/h_norm * e_norm * sin(theta_star);
v_theta = mu_mars/h_norm * (1 + e_norm*cos(theta_star));
v2 = [v_r v_theta 0];
v2_norm = norm(v2);
v2_hat = v2/v2_norm;

% Part c
theta_star_p = 0;
theta_star_a = pi;

theta_p = theta_star_p + aop;
C_XYZ_RTH_p = R3(-raan) * R1(-i) * R3(-theta_p);

r_p = a * (1 - e_norm);
r_p_rth = [r_p 0 0];
v_r_p = mu_mars/h_norm * e_norm*sin(theta_star_p);
v_theta_p = mu_mars/h_norm * (1 + e_norm*cos(theta_star_p));
v_p_rth = [v_r_p v_theta_p 0];
r_p_xyz = C_XYZ_RTH_p * r_p_rth';
v_p_xyz = C_XYZ_RTH_p * v_p_rth';

r_a = a * (1 + e_norm);
r_a_rth = [r_a 0 0];
v_r_a = mu_mars/h_norm * e_norm*sin(theta_star_a);
v_theta_a = mu_mars/h_norm * (1 + e_norm*cos(theta_star_a));
v_a_rth = [v_r_a v_theta_a 0];
r_a_xyz = C_XYZ_RTH * r_a_rth';
v_a_xyz = C_XYZ_RTH * v_a_rth';











