close all; clear; clc;

% ASEN 6020, HW 1, Problem 1
% Jash Bhalavat
% Spring 2025

n = 10000;
r_range_gt1 = linspace(1, 15, n);
r_range_lt1 = linspace(0, 1, n);

for i = 1:n
    cost_gt1(i) = jbp_vs_jh(r_range_gt1(i));
    cost_lt1(i) = r_less_than_1(r_range_lt1(i));
    cost(i) = cost_compare(r_range_lt1(i));
end

figure()
plot(r_range_gt1, cost_gt1)

figure()
plot(r_range_lt1, cost_lt1)

function out = jbp_vs_jh(r)
    out = sqrt(r) - (sqrt(2) - 1) - (r-1)/(sqrt(1+r));
end

function out = cost_compare(r)
    jbp = (sqrt(2)-1) * (1 + 1/sqrt(r));
    jh = sqrt((2*r)/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)));
    out = jbp < jh;
end

function out = r_less_than_1(r)
    a = sqrt(2) - 1;
    out = 5*r^2 + 6*r + (1+r)*(4*sqrt(r)*sqrt(1+r) - (2*a + a*sqrt(1+r))^2);
end