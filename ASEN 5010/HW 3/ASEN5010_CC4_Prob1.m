clear; clc; close all;

% ASEN 5010 - HW 3, Concept Check 4 KE, Problem 1
% Spring 2025
% Jash Bhalavat

m1 = 1;
m2 = 1;
m3 = 2;
m4 = 2;
tot_mass = m1 + m2 + m3 + m4;

R1 = [1 -1 2];
R2 = [-1 -3 2];
R3 = [2 -1 -1];
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

H_o = m1*cross(R1, R1d) + m2*cross(R2, R2d) + m3*cross(R3, R3d) + m4*cross(R4, R4d);

H_c = m1*cross(r1, r1d) + m2*cross(r2, r2d) + m3*cross(r3, r3d) + m4*cross(r4, r4d);

% H0 = [0, -4, 18] # Nms
    % Hc = [1.33333333333333, 2, 0.666666666666667] # Nms 
