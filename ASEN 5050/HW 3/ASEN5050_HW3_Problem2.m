clear; clc; close all;

% GMAT
mu_mars = 42828.314258067;

a = 6463.8;
e = 0.45454;
i = deg2rad(74.924);
raan = deg2rad(1.241);
aop = deg2rad(353.31);
theta_star = deg2rad(199.38);

period = 2*pi * sqrt(a^3/mu_mars);
r_p = a * (1 - e);

% Part c
coe = [a, e, i, raan, aop - 2*pi, theta_star - 2*pi];

rv = coe2rv(coe, mu_mars);

% Part g
ecc_point_mass = readmatrix("Problem3_ECC");
hmag_point_mass = readmatrix("Problem3_HMAG");

figure()
plot(ecc_point_mass(:,1), ecc_point_mass(:,2))
title("Mars Point Mass Model Eccentricity")
xlabel("Days")
ylabel("Eccentricity")

figure()
plot(hmag_point_mass(:,1), hmag_point_mass(:,2))
title("Mars Point Mass Model Specific Angular Momentum")
xlabel("Days")
ylabel("Specific Angular Momentum")

% Part j
ecc_higher_fidelity = readmatrix("Problem3_ECC_46");
hmag_higher_fidelity = readmatrix("Problem3_HMAG_46");

figure()
plot(ecc_higher_fidelity(:,1), ecc_higher_fidelity(:,2))
title("Mars Higher Fidelity Model Eccentricity")
xlabel("Days")
ylabel("Eccentricity")

figure()
plot(hmag_higher_fidelity(:,1), hmag_higher_fidelity(:,2))
title("Mars Higher Fidelity Model Specific Angular Momentum")
xlabel("Days")
ylabel("Specific Angular Momentum")

