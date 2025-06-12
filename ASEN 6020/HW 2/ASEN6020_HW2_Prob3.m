clear; clc; close all;

% ASEN 6020 - HW 2, Problem 3
% Spring 2025
% Jash Bhalavat

n = 100;

function out = f(x,y)
    out = -sin(x)*cos(y);
end

xrange = linspace(-10, 10, n);
yrange = linspace(-10, 10, n);

for i = 1:n
    for j = 1:n
        z(i,j) = f(xrange(i), yrange(j));
    end
end

figure()
contourf(xrange, yrange, z')
xlim([0, pi])
ylim([0, pi])
colorbar
