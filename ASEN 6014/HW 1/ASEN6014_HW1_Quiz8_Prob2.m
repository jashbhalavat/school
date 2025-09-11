clear; clc;

% ASEN 6014 - HW 1, Quiz 8, Prob 2
% Fall 2025
% Jash Bhalavat

a = 8000.0; % km
v = [-6.80822,1.04998,3.61939]; % km/s
v_sq = norm(v)^2;
mu = 398600;

% vis-viva eqn.
r = 2 / (v_sq/mu + 1/a); % km