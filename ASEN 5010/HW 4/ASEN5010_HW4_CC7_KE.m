clear; clc; close all;

I_c = [10 1 -1;
        1 5 1;
        -1 1 8];
omega_B = [0.01; -0.01; 0.01];

T = 1/2 * omega_B' * I_c * omega_B;