clear; clc; close all;

% ASEN 5050, HW 5, Problem 3
% Fall 2024, 10/17/2024
% Jash Bhalavat

mu_sun = 1.32712428e11;
mu_saturn = 3.794e7;

AU = 149597870.7;

a_earth = 1.0000010178 * AU;
a_saturn = 9.554909595 * AU;

% Given

% Part a
% Earth to saturn hohmann transfer
hohmann_out = hohmann_transfer(a_earth, a_saturn, mu_sun);
delta_v_hohmann = hohmann_out(1);
tof_hohmann = hohmann_out(2);


% Part b
n2 = sqrt(mu_sun / (a_saturn)^3);
alpha = n2 * tof_hohmann;

% Phasing angle
phi = pi - alpha;


% Par c, Bi-elliptic transfer
r_b = 11 * AU;

bi_elliptic_out = bi_elliptic_transfer(a_earth, a_saturn, r_b, mu_sun);
delta_v_be = bi_elliptic_out(1);
tof_be = bi_elliptic_out(2);


