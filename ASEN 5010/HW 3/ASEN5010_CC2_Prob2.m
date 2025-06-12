clear; clc; close all;

% ASEN 5010 - HW 3, Concept Check 1 KE, Problem 2
% Spring 2025
% Jash Bhalavat

m1 = 1;
m2 = 1;
m3 = 2;
m4 = 2;
tot_mass = m1 + m2 + m3 + m4;

R1 = [1 -1 2];
R2 = [-1 -3 2];
R3 = [2 -1 1];
R4 = [3 -1 -2];

Rc = (m1*R1 + m2*R2 + m3*R3 + m4*R4)/tot_mass;
r1 = R1 - Rc;
r2 = R2 - Rc;
r3 = R3 - Rc;
r4 = R4 - Rc;

R1d = [2 1 1];
R2d = [0 -1 1];
R3d = [3 2 -1];
R4d = [0 0 1];

Rcd = (m1*R1d + m2*R2d + m3*R3d + m4*R4d)/tot_mass;
r1d = R1d - Rcd;
r2d = R2d - Rcd;
r3d = R3d - Rcd;
r4d = R4d - Rcd;

transEnergy = 1/2 * tot_mass * dot(Rcd, Rcd);
rotEnergy = 1/2 * (dot(r1d, r1d)*m1 + dot(r2d, r2d)*m2 + dot(r3d, r3d)*m3 + dot(r4d, r4d)*m4);

% transEnergy = 7 # Joules
    % rotEnergy = 12 #Joules
