clear; clc;

mu_sun = 1.32712428e11;

AU = 149597870.7;

mars_orbit_radius = 1.52*AU;

a = 5.75*AU;
e = 0.8104;
p = a * (1 - e^2);
h = sqrt(p * mu_sun);
b = a * sqrt(1-e^2);

r_intersection_w_titan = mars_orbit_radius;
theta_star = acos(1/e * (p/r_intersection_w_titan - 1));

v_r = mu_sun/h * e * sin(theta_star);
v_theta = mu_sun/h * (1 + e * cos(theta_star));
v_sc_intersection = [v_r, v_theta, 0];

v_titan_theta = sqrt(mu_sun / mars_orbit_radius);
v_titan = [0 v_titan_theta 0];
v_titan_mag = norm(v_titan);

v_inf_in = v_sc_intersection - v_titan;