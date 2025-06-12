close all; clear; clc;

% ASEN 5010 - HW 2, Problem 3 - 3.17
% Spring 2025
% Jash Bhalavat

syms phi a b c

function out = skew(x)
    out = [0 -x(3) x(2);
            x(3) 0 -x(1);
            -x(2) x(1) 0];
end

% test = expm(-phi*skew([a b c]));
x = [a b c];
x_squared = a^2 + b^2 + c^2;
(skew(x) * skew(x) * skew(x) * skew(x))

