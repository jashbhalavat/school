clear; clc; close all;

syms t;

A = [0 1 0; 0 0 1; 0 0 0]

[l, v] = eig(A)

eat = expm(A*t)