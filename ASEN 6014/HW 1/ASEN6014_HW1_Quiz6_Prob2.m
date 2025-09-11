clear; clc; close all;

% ASEN 6014 - HW 1, Quiz 6, Problem 2, 3
% Fall 2025
% Jash Bhalavat

r = [2466.69,5941.54,3282.71]; % km
v = [-6.80822,1.04998,3.61939]; % km/s
h = cross(r, v); 

h_mag = norm(h);
radius = 8000; % km
f_dot = (h_mag/radius^2) * 180/pi; % deg/s

% # adjust the return matrix values as needed
% def result():
%     hVec_N = [18057.9706148, -31277.3249953, 43041.286625]  # km
% 
%     return hVec_N




