clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 3, Problem 3
% Jash Bhalavat

%% Given

x = 10; % km
y = 500; % km
r = 7000; % km
xdot = 0.1; % km/s
ydot = -0.1; % km/s
rdot = 0.05; % km/s

[delta_r, delta_r_dot, s, s_dot] = rect2curv(x, y, r, xdot, ydot, rdot);

% # adjust the return matrix values as needed
% def result():
%     delta_r = [27.8090469220915]  # km
%     delta_r_dot = [0.092505294795768] # km/s
%     s = [498.442602242367] # km
%     s_dot = [-0.106421250706256] # km/s
%     return delta_r, delta_r_dot, s, s_dot

function [delta_r, delta_r_dot, s, s_dot] = rect2curv(x, y, r, xdot, ydot, rdot)
    rd = sqrt((x+r)^2+y^2);
    delta_r = rd - r;
    s = r * atan(y/(r+x));
    % delta_r_dot = ((r+x)*(rdot+xdot + y*ydot))/rd;
    rd_dot = 1/2 * ((r+x)^2 + (y)^2)^(-1/2) * (2*(r+x)*(rdot + xdot) + 2*y*ydot);
    delta_r_dot = rd_dot - rdot;
    % s_dot = ((ydot*r - y*rdot)/r)/(1 + (y/r)^2) + (atan(y/r)*rdot);
    % s_dot = (ydot*r - y*rdot)/(r + y^2/r) + atan(y/r)*rdot;
    % s_dot = r * ((ydot*r - y*rdot)/(r^2))/(1 + (y/r)^2) + atan(y/r)*rdot;
    delta_theta_dot = (1/(1+(y/(r+x))^2)) * ((ydot*(r+x) - y*(rdot + xdot))/(r+x)^2);
    delta_theta = atan((y)/(r+x));
    s_dot = r*delta_theta_dot + rdot*delta_theta;

end
