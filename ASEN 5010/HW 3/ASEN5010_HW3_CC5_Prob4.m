clear; clc; close all;

% ASEN 5010 - CC5 Rigid Body Kinematics, Problem 4
% Spring 2025
% Jash Bhalavat

D2R = pi/180;
t1 = -10*D2R;
t2 = 10*D2R;
t3 = 5*D2R;

BN = R1(t3)*R2(t2)*R3(t1);

omega_BN_N = [0.01, -0.01, 0.01];

I_B = [10, 1, -1;
    1, 5, 1;
    -1, 1, 8];

omega_BN_B = BN * omega_BN_N';
H = I_B * omega_BN_B;
% H = [0.077152176539363, -0.013041792460219, 0.083453291546846] # Nms 
