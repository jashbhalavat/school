clear; clc;

% ASEN 6014 - HW 1, Quiz 10, Prob 6
% Fall 2025
% Jash Bhalavat

mu = 398600;
a = 8000; % km
e = 0.1;
i = 30.0 * pi/180; % rad
omega = 120.0 * pi/180; % rad
Omega = 145 * pi/180; % rad
M_t0 = 10 * pi/180.0; % rad
P = 2*pi*sqrt(a^3/mu); % sec

t0 = M_t0/(2*pi)*P; % sec
t1 = t0 + 3600; % sec

M_t1 = (t1*2*pi)/P; % rad

% Root solving
fun = @(E) E - e*sin(E) - M_t1;
E_t1 = fzero(fun, M_t1);

f_t1 = 2 * atan(sqrt((1+e)/(1-e)) * tan(E_t1/2)); % rad

rv_t1 = coe2rv([a, e, i, Omega, omega, f_t1], mu);

% # adjust the return matrix values as needed
% def result():
%     rVec = [-1264.60766618468, 8013.80931615287, -3371.25163984785]  # km
%     vVec = [-6.03962093013619, -0.204397604896008, 2.09671503286582]  # km/sec
% 
%     return rVec, vVec
