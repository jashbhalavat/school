clear; clc; close all;

T = [1 0 -1 0; 0 1 0 -1; 1 0 1 0; 0 1 0 1];
T_inv = inv(T);

syms k1 k2 k3 m1 m2;

A = [0 1 0 0; (-k1-k2)/m1 0 k2/m1 0; 0 0 0 1; k2/m2 0 (-k2-k3)/m2 0];
B = [0 0; -1/m1 0; 0 0; 1/m1 1/m2];

A_tilde = T * A * T_inv;
B_tilde = T * B;



