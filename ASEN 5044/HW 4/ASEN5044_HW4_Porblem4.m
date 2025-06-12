clear; clc; close all;

% Problem 4
% Textbook - 2.16
n = 10000;

lower_bound = -1/2;
upper_bound = 1/2;

x1 = lower_bound + (upper_bound - lower_bound).*rand(n, 1);
x2 = lower_bound + (upper_bound - lower_bound).*rand(n, 1);
x3 = lower_bound + (upper_bound - lower_bound).*rand(n, 1);
x4 = lower_bound + (upper_bound - lower_bound).*rand(n, 1);

samples12 = (x1 + x2)/2;
samples1234 = (x1 + x2 + x3 + x4)/4;

nbins = 50;

figure()
histogram(samples12, nbins)
xlabel("(x1+x2)/2")
ylabel("Frequency")
title("Histogram for (x1 + x2)/2")

figure()
histogram(samples1234, nbins)
xlabel("(x1+x2+x3+x4)/4")
ylabel("Frequency")
title("Histogram for (x1 + x2 + x3 + x4)/4")



