clear; clc; close all;

% ASEN 6060 - HW 1, Problem 1
% Jash Bhalavat

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

%% Part a
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






