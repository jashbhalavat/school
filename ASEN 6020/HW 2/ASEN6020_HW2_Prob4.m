clear; clc; close all;

% ASEN 6020 - HW 2, Problem 4
% Spring 2025
% Jash Bhalavat

n = 100;

function out = f(x,y)
    x_vec = [x; y];
    out = x_vec' * [1 1; 1 1] * x_vec;
end

xrange = linspace(-1, 1, n);
yrange = linspace(-1, 1, n);

for i = 1:n
    for j = 1:n
        z(i,j) = f(xrange(i), yrange(j));
    end
end

figure()
contourf(xrange, yrange, z')
xlabel("x1")
ylabel("x2")
title("No Unique Minimum when Q is semi-definite")
colorbar
