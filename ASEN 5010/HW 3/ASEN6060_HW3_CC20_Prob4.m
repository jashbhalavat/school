clear; clc; close all;

% ASEN 6060 - HW 3, CC20, Problem 4
% Spring 2025
% Jash Bhalavat

sigma0 = [0.4, 0.2, -0.1];

tspan = [0:60]; % sec

% [t, sigma] = ode45(@diff_eq, [0:60], sigma0);
% [t, sigma] = rkf45(@diff_eq, tspan, sigma0);
sigma = rk4(@diff_eq, [0:60], sigma0);

norm_at_42 = sqrt(sigma(43, 1)^2 + sigma(43, 2)^2 + sigma(43, 3)^2);

function out = omega_B(t)
    out = 20 * [sin(0.1*t); 0.01; cos(0.1*t)] * pi/180;
end

function out = B(mrp)
    mrp_sq = norm(mrp)^2;
    out = (1-mrp_sq)*eye(3) + 2*skew_symmetric(mrp) + 2*mrp*mrp';
end

function ydot = diff_eq(t, y)
    if norm(y) > 1
        y = -y./norm(y)^2;
    end
    omega_B_at_t = omega_B(t);
    ydot = 1/4*B(y) * omega_B_at_t;
end

function out = skew_symmetric(vec)
    out = [0 -vec(3) vec(2);
            vec(3) 0 -vec(1);
            -vec(2) vec(1) 0];
end

