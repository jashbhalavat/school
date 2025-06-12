close all; clear; clc;

% ASEN 6020, HW 1, Problem 2
% Jash Bhalavat
% Spring 2025

r_range = linspace(1,500,500);

cost = zeros(500) + 2;

j_be(10, 5)
j_h(10)

for i = 1:length(r_range)
    l_range = r_range(1:i);
    for j = 1:length(l_range)
        cost(i,j) = j_be(r_range(i), l_range(j)) < j_h(r_range(i));
    end
end

figure()
contourf(r_range, l_range, cost')
colorbar();
title("Cost of Bi-Elliptic vs Hohmann Transfer")
xlabel("r", 'FontSize',15)
ylabel("l", 'FontSize',15)

function out = j_be(r, l)
    out = sqrt((2*l)/(1+l)) - 1 + sqrt((2*r)/(l*(l+r))) - sqrt(2/(l*(l+1))) + 1/sqrt(r) - sqrt((2*l)/(r*(l+r)));
end

function out = j_h(r)
    out = sqrt((2*r)/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)));
end
