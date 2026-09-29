clear; clc; close all;

%% ASEN 5090 - HW 2
% Jash Bhalavat
% 09/12/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

% Earth gravitational constant
mu = 3.986005e14; % m3/s2, Lecture 5, slide 31

% Earth rotation rate
omega_dot = 7.2921151467e-5; % rad/s, Lecture 5, slide 31

%% Add path

addpath("HW2_Code/")
addpath("HW2_DATA/")

%% Problem 5c

time = (0:100:604800)';

filename = 'YUMA231.alm.txt';
ephem_all = read_GPSyuma(filename);
week_number = 384 * ones([length(time), 1]);
t_input =  [week_number, time];
[prn05_alm_health, prn05_alm_pos, ~] = alm2pos(ephem_all, t_input, 5);

figure(1);
subplot(3,1,1)
plot(time, prn05_alm_pos(:,1), '--', 'LineWidth', 2)
xlim([0, time(end)])
grid on
ylabel("X [m]")

subplot(3,1,2)
plot(time, prn05_alm_pos(:,2), '--', 'LineWidth', 2)
grid on
xlim([0, time(end)])
ylabel("Y [m]")

subplot(3,1,3)
plot(time, prn05_alm_pos(:,3), '--', 'LineWidth', 2)
grid on
xlim([0, time(end)])
ylabel("Z [m]")
xlabel("Time of Week 2432 [s]")
sgtitle("PRN05 ECEF Positions predicted from YUMA almanac")

%% Remote path

rmpath("HW2_Code/")
rmpath("HW2_DATA/")