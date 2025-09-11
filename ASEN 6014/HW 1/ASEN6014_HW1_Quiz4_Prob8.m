clear; clc; close all;

% ASEN 6014 - HW 1, Quiz 4, Problem 8
% Fall 2025
% Jash Bhalavat

a = 8000.0; % km
e = 0.23;
f = 295.0 * pi/180.0; % rad

rp = a * (1 - e); % km
ra = a * (1 + e); % km
p = a * (1 - e^2); % km
b = a * sqrt(1 - e^2); % km
E = ((2 * atan(sqrt((1-e)/(1+e)) * tan(f/2))) + 2*pi) * 180.0/pi; % deg
