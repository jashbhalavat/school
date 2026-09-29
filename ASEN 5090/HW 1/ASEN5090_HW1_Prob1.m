clear; clc; close all;

%% ASEN 5090 - HW 1
% Jash Bhalavat
% 08/24/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

%% Problem 1a

rho_PL1 = 550; % m
rho_PL2 = 500; % m

delta_r = (rho_PL1 - rho_PL2) / (2*c); % sec
R_r_PL1 = rho_PL1 - c*delta_r; % m
R_r_PL2 = 1000 - R_r_PL1; % m

