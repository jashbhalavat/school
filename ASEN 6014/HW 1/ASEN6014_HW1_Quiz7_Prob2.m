clear; clc;

% ASEN 6014 - HW 1, Quiz 7, Prob 2
% Fall 2025
% Jash Bhalavat

r_t0 = [2466.69,5941.54,3282.71]; % km
v_t0 = [-6.80822,1.04998,3.61939]; % km/s
v_t = [5.57433,-0.92203,-3.00873]; % km/s

mu = 398600; % km3/s2

h = cross(r_t0, v_t0); % constant
p = norm(h)^2/mu;
e = (cross(v_t0, h))/(mu) - r_t0/norm(r_t0); % e = (v x h)/mu - r/|r|

e_unit = e/norm(e); % constant
rt_unit = (cross(v_t, h))/(mu) - e; % r/|r| = r_hat = e - (v x h)/mu
f = acos(dot(e_unit, rt_unit)); % f = ang diff --> e, r
rt_norm = p/(1+norm(e)*cos(f)); % r = p/(1+ecos(f))
rt = rt_unit * rt_norm;

% # adjust the return matrix values as needed
% def result():
%     rHat_N = [-0.391500240332171, -0.814701127362729, -0.427774668061654]
% 
%     return rHat_N



