clear; clc; close all;

%% Call varying mass function

n = 25;

final_nondim_mass = linspace(0.95, 0.96, n);

for i = 2:length(final_nondim_mass)
    output(i,:) = varying_final_mass(final_nondim_mass(i));
end

