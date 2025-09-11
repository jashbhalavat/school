clear; clc;

% ASEN 6014 - HW 1, Quiz 9, Prob 3
% Fall 2025
% Jash Bhalavat

mu = 398600;

a = 7500; % km
P = 2*pi*sqrt(a^3/mu); % sec
n = (2*pi)/P; % rad/s
e = 0.05;
f_t0 = 25*pi/180; % rard
E_t0 = 2 * atan(sqrt((1-e)/(1+e)) * tan(f_t0/2)); % rad
t0 = abs((E_t0 - e*sin(E_t0))/n); % sec

tf = t0 + 3600; % s after t0, sec
E_tf = kepler_solver_eclipse(tf, a, e, mu); % rad
f_tf = 2 * atan(sqrt((1+e)/(1-e)) * tan(E_tf/2)); % rad
if f_tf < 0
    f_tf = f_tf + 2*pi;
end



