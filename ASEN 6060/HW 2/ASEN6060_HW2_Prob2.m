clear; clc; close all;

% ASEN 6060 - HW 2, Problem 2
% Jash Bhalavat

%% Constants

G = 6.67408 * 10^-11; % m3/(kgs2)
G = G / (10^9); % km3/(kgs2)

% Earth
mu_earth = 398600.435507; % km3/s2
a_earth = 149598023; % km
e_earth = 0.016708617;
mass_earth = mu_earth / G; % kg

% Moon
mu_moon = 4902.800118; % km3/s2
a_moon = 384400; % km
e_moon = 0.05490;
mass_moon = mu_moon / G; % kg

% Sun
mu_sun = 132712440041.279419; % km3/s2
mass_sun = mu_sun / G; % kg

% Earth-Moon system
mass_ratio_em = mass_moon / (mass_earth + mass_moon);
m_star_em = mass_earth + mass_moon;
l_star_em = a_moon;
t_star_em = sqrt(l_star_em^3/(G * m_star_em));

% Sun-Earth system
mass_ratio_se = mass_earth / (mass_earth + mass_sun);
m_star_se = mass_earth + mass_sun;
l_star_se = a_earth;
t_star_se = sqrt(l_star_se^3/(G * m_star_se));

%% Evals for EM system

mu = mass_ratio_em;

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Calculate out of plane modes for all 5 eq points
for i = 1:5
    em_eq_pts_out_of_plane_modes(i,:) = out_of_plane_modes(mu, em_eq_pts(i,:));
    em_eq_pts_in_plane_modes(i,:) = in_plane_modes(mu, em_eq_pts(i,:));
end

%% Evals for SE system

mu = mass_ratio_se;

% Earth Moon system equilibrium points
[se_eq_pts, se_eq_validity] = all_eq_points(mu);

% Calculate out of plane modes for all 5 eq points
for i = 1:5
    se_eq_pts_out_of_plane_modes(i,:) = out_of_plane_modes(mu, se_eq_pts(i,:));
    se_eq_pts_in_plane_modes(i,:) = in_plane_modes(mu, se_eq_pts(i,:));
end

