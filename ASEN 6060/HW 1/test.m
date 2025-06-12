clear; clc; close all;

% ASEN 6060 - HW 1, Problem 4
% Jash Bhalavat
% 01/28/2025

%% Constants
G = 6.67408 * 10^-11; % m3/(kgs2)
G = G / (10^9); % km3/(kgs2)
% Earth
mu_earth = 398600.435507; % km3/s2
mass_earth = mu_earth / G; % kg
% Moon
mu_moon = 4902.800118; % km3/s2
a_moon = 384400; % km
mass_moon = mu_moon / G; % kg
% Earth-Moon system
mass_ratio_em = mass_moon / (mass_earth + mass_moon);
m_star_em = mass_earth + mass_moon;
mu = mass_ratio_em;

%% Zero Velocity Surfaces

mu = mass_ratio_em;
n = 5001;
x = linspace(-2, 2, n);
y = linspace(-2, 2, n);

c_given = 3.189;

zvc = full_zvc(c_given, mu, x, y);

function zvc = full_zvc(c_given, mu, x, y)
    % Initialize zero-velocity curve as a ones matrix
    zvc = ones([length(x), length(y)]);
    
    % Parse through x and y and calculate 2U_star at each point
    for i = 1:length(x)
        for j = 1:length(y)
            c_calc = u_star_times_2(x(i), y(j), mu);
            % If 2U_star is less than C, then that point is not allowable
            % motion. Assign zvc at that point to -1
            if c_calc < c_given
                zvc(i, j) = -1;
            end
        end
    end
end

function out = u_star_times_2(x, y, mu)
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2);
    out = (x^2 + y^2) + 2*(1 - mu)/r1 + 2*mu/r2;
end


figure()
contourf(x, y, zvc')
map = [220/255, 220/255, 220/255
    1, 1, 1];
colormap(map)
hold on
scatter(-mu, 0, 200, 'filled')
scatter(1-mu, 0, 50, 'filled')
hold off
legend("", "Earth", "Moon")
xlabel("x [Non-Dimensional]")
ylabel("y [Non-Dimensional]")
title("Zero velocity curve for Jacobi Constant C = " + num2str(c_given))