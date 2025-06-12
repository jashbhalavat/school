clear; clc; close all;

% Given
mu_saturn = 3.794e7;
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
theta_star_1 = coe(6);

h = cross(r1, v1);
h_norm = norm(h);

% Part a
% Point 2 is impact on surface
theta_star_2 = abs(acos(((h_norm^2/mu_saturn)/saturn_r_eq - 1)/e)) * -1; 

r2v2 = fg(r1, v1, theta_star_2, mu_saturn);

% Part b
P = 2 * pi * sqrt(a^3/mu_saturn);
n = 2*pi/P;

E_1 = atan(sqrt((1-e)/(1+e))*tan(theta_star_1/2)) * 2;
E_2 = atan(sqrt((1-e)/(1+e))*tan(theta_star_2/2)) * 2;

if E_1 < 0
    E_1 = E_1 + 2*pi;
end

if E_2 < 0
    E_2 = E_2 + 2*pi;
end

t1 = (E_1 - e*sin(E_1))/n;
t2 = (E_2 - e*sin(E_2))/n;

time_taken = t2 - t1;







