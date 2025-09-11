clear; clc;

% ASEN 6014 - HW 1, Quiz 10, Prob 7
% Fall 2025
% Jash Bhalavat

mu = 398600;

r_N = [-820.865,-1905.95,-7445.9];
v_N = [-6.75764,-1.85916,0.930651];

coe = rv2coe(r_N, v_N, mu);

% # adjust the return matrix values as needed
% def result():
%     sma = [7499.99576967534]  # km
%     ecc = [0.0500000171539765]  
%     inc = [1.78023546818388]  # rad
%     AN = [-2.84488718870787]   # rad
%     AP = [2.61798287192296]   # rad
%     f = [2.26893877004198]    # rad
% 
%     return sma, ecc, inc, AN, AP, f
