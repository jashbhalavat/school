clear; clc; close all;

% Given
mu_jupiter = 1.268e8;
jupiter_r_eq = 71492;

% t1 - Aug 29, 1996
r1 = [5.352950e6, 7.053778e5, -4.059700e5];
v1 = [-4.164248, 1.963690, 3.191257e-1];

% Part a
coe = rv2coe(r1, v1, mu_jupiter);
theta_star_1 = coe(6);
e = coe(2);
E_1 = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta_star_1/2));

% Part b
% t2 - s/c crosses the descending node
aop = coe(5);
theta_star_2 = pi - aop;
E_2 = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta_star_2/2));

% Part c
a = coe(1);
P = 2*pi*sqrt(a^3/mu_jupiter);
n = (2*pi)/P;

% since E and theta_star are negative, this time is time until periapsis
t1 = abs((E_1 - e*sin(E_1))/n);

% since E and theta_star are positive, this time is time from periapsis
t2 = (E_2 - e*sin(E_2))/n;

% t1 to t2 in days
t1_to_t2 = (t1 + t2)/(24*60*60);

% Part d
r2v2 = fg(r1, v1, theta_star_2, mu_jupiter);

% Part f
% t3 is 20 days after t1
% t1 indicates time until periapsis
% t3 time since periapsis is 20 - t1
t3_time_since_periapsis = 20*86400 - t1;

E_3 = kepler_solver_eclipse(t3_time_since_periapsis, a, e, mu_jupiter);

theta_star_3 = 2 * atan(sqrt((1+e)/(1-e)) * tan(E_3/2));


