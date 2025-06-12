clear; clc; close all;

% ASEN 5044 - Final Project - Stat OD
% Jash Bhalavat

% Data
data = load('orbitdeterm_finalproj_KFdata.mat');

% Constants
mu = 398600; % km^3/s^2
r0 = 6678; % km
x0 = r0; % km
y0 = 0; % km
xdot0 = 0; % km/s
ydot0 = r0 * sqrt(mu/r0^3); % km/s
n = 4; % state variables
m = 2; % control variables
R_E = 6378; % km, earth radius
omega_E = (2*pi)/86400; % rad/s
h = 300; % km, satellite altitude

% Satellite orbit period [sec]
P = 2*pi*sqrt((h + R_E)^3/mu);

% Let T = Orbit period [sec]
T = P;

% Ground station initial location
number_of_ground_stations = 12;
for i = 1:number_of_ground_stations
    theta_i_0(i, 1) = (i-1)*pi/6;
    X_i_0(i, 1) = R_E * cos(theta_i_0(i));
    Y_i_0(i, 1) = R_E * sin(theta_i_0(i));
end

% Store in array
theta_i = [theta_i_0];
X_i = [X_i_0];
Y_i = [Y_i_0];

% Calculate for time T. Make sure that time is in time units and not
% discrete time steps
for time = 1:T
    for i = 1:number_of_ground_stations
        X_i(i, time+1) = R_E * cos(theta_i_0(i) + omega_E*time);
        Y_i(i, time+1) = R_E * sin(theta_i_0(i) + omega_E*time);
        theta_i(i, time+1) = atan(Y_i(i, time+1)/X_i(i, time+1));
    end
end

ground_station_elevation_range_lower = theta_i - pi/2;
ground_station_elevation_range_upper = theta_i + pi/2;

time_vec = 0:T;
figure()
plot(time_vec, X_i(1,:), time_vec, Y_i(1,:))

