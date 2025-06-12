clear; clc; close all;

beta0 = [0.408248,0.,0.408248,0.816497]';

tspan = [0:60]; % sec

[t, beta] = ode45(@diff_eq, tspan, beta0);

norm_at_42 = sqrt(beta(43, 2)^2 + beta(43, 3)^2 + beta(43, 4)^2);

function out = omega_B(t)
    out = 20 * [sin(0.1*t); 0.01; cos(0.1*t)] * pi/180;
end

function out = mult_mat(omega_in)
    out = [0, -omega_in(1), -omega_in(2), -omega_in(3);
            omega_in(1), 0, omega_in(3), -omega_in(2);
            omega_in(2), -omega_in(3), 0, omega_in(1);
            omega_in(3), omega_in(2), -omega_in(1), 0];
end

function ydot = diff_eq(t, y)
    omega_B_at_t = omega_B(t);
    ydot = 1/2*mult_mat(omega_B_at_t) * y;
end

