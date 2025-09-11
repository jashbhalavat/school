clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 3, Problem 4
% Jash Bhalavat

%% Given

delta_r = 10; % km
s = 500; % km
r = 7000; % km
delta_r_dot = 0.1; % km/s
s_dot = -0.1; % km/s
rdot = 0.05; % km/s

[x, y, xdot, ydot] = curv2rect(delta_r, s, r, delta_r_dot, s_dot, rdot);

% # adjust the return matrix values as needed
% def result():
%     x = [-7.91538989652236]  # km
%     x_dot = [0.107585010328153] # km/s
%     y = [500.852079115282] # km
%     y_dot = [-0.100524141953074] # km/s
%     return x, x_dot, y, y_dot

function [x, y, xdot, ydot] = curv2rect(delta_r, s, r, delta_r_dot, s_dot, rdot)
    % y = (r+x) * tan(s/r);
    % x = sqrt((delta_r + r)^2 - y^2) - r;
    % ydot = (r+x)*(sec(s/r))^2*((r*s_dot - s*rdot)/(r^2)) + tan(s/r)*(rdot+xdot);
    % xdot = 1/2*((delta_r + r)^2 - y^2)^(-1/2)*(2*(delta_r + r)*(delta_r_dot + rdot) - 2*y*ydot) - rdot;

    sec_sr_sq = (sec(s/r))^2;
    x = (-2*r*sec_sr_sq + sqrt(4*r^2*sec_sr_sq^2 - 4*sec_sr_sq*(r^2*(tan(s/r))^2 - delta_r^2 - 2*r*delta_r)))/(2*sec_sr_sq);
    y = (r+x)*tan(s/r);

    a = (r*s_dot - s*rdot)/(r^2);

    num_der_1 = -2*r*2*sec_sr_sq*tan(s/r)*a - sec_sr_sq*2*rdot;
    num_der_2 = 1/2 * (4*sec_sr_sq*(-delta_r^2 - 2*r*delta_r))^(-1/2) * ((4*sec_sr_sq)*(-2*delta_r_dot - (2*r*delta_r_dot + 2*delta_r*rdot)) + (-delta_r^2 - 2*r*delta_r)*(8*sec_sr_sq*tan(s/r)*a));

    xdot = ((2*sec_sr_sq)*(num_der_1 + num_der_2) - (-2*r*sec_sr_sq + sqrt(4*r^2*sec_sr_sq^2 - 4*sec_sr_sq*(r^2*(tan(s/r))^2 - delta_r^2 - 2*r*delta_r)))*(4*sec_sr_sq*tan(s/r)*a))/(4*sec_sr_sq);

end