close all; clear; clc;

% ASEN 6020, HW 1, Problem 6
% Jash Bhalavat
% Spring 2025

e_range = linspace(0, 1, 1000);

for i = 1:length(e_range)
    dva(i) = e_range(i);
    dvb(i) = sqrt(1-e_range(i)) - 1 + e_range(i);
    dvc(i) = 1 + e_range(i) - sqrt(1 + e_range(i));
end

figure()
plot(e_range, dva, 'LineWidth',2.0)
hold on
plot(e_range, dvb, 'LineWidth',2.0)
plot(e_range, dvc, 'LineWidth',2.0)
hold off
legend("One-Impulse", "Two-Impulse Apoapsis", "Two-Impulse Periapsis")
title("Normalized Cost for a 180° AOP change")
xlabel("e")
ylabel("Normalized Cost")
