close all; clear; clc;

% ASEN 5010 - HW 2, Problem 1 - 3.13
% Spring 2025
% Jash Bhalavat

D2R = pi/180;

t1 = -30*D2R;
t2 = 40*D2R;
t3 = 20*D2R;

C = R3(t3) * R1(t2) * R3(t1);

%% Part a, b
phi = acos(1/2 * (C(1,1) + C(2,2) + C(3,3) - 1));
e_hat = 1/(2*sin(phi)) * [C(2,3)-C(3,2), C(3,1)-C(1,3), C(1,2)-C(2,1)];
phi_prime = 2*pi - phi;

%% Part c
ep = sheppards_method(C);

%% Part d
crp = [ep(2)/ep(1), ep(3)/ep(1), ep(4)/ep(1)];
crp2 = tan(phi/2) * e_hat;

%% Part e
mrp = [ep(2)/(1+ep(1)), ep(3)/(1+ep(1)), ep(4)/(1+ep(1))];
mrp2 = tan(phi/4) * e_hat;


