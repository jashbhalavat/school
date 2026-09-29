clear; clc; close all;

%% ASEN 5090 - HW 1
% Jash Bhalavat
% 08/24/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

n = 500;

%% Problem 2

v0 = 50; % m/s
h0 = 100; % m
x0 = 250; % m

xA = linspace(-x0, x0, n);

for i = 1:n
    rho(i) = sqrt(h0^2 + xA(i)^2); % m
    rho_dot(i) = v0*xA(i) / rho(i); % m/s
    z(i) = atan2(xA(i), h0);
end

% Plot
figure(1)
subplot(3,1,1)
plot(xA, rho, 'LineWidth', 2)
ylabel('$\rho$ [m]', 'Interpreter', 'Latex', 'FontSize', 15)
grid on

subplot(3,1,2)
plot(xA, rho_dot, 'LineWidth', 2)
ylabel('$\dot{\rho}$ [m/s]', 'Interpreter', 'Latex', 'FontSize', 15)
grid on

subplot(3,1,3)
plot(xA, z, 'LineWidth', 2)
ylabel('$z$ [rad]', 'Interpreter', 'Latex', 'FontSize', 15)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15)
grid on

sgtitle("Range, Range Rate, and Zenith Angle Measurements")



