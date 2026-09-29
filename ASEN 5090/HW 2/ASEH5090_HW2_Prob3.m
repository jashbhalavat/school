clear; clc; close all;

%% ASEN 5090 - HW 2
% Jash Bhalavat
% 09/11/2026

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

%% Problem 3a

% Obtain SP3 data for all sats
filename = 'IGS0OPSFIN_20262310000_01D_15M_ORB.SP3';
sp3_data = read_sp3(filename);

% 
prn04_sp3 = sp3_data(sp3_data(:, 3) == 4, :);
time = (prn04_sp3(:,2) - prn04_sp3(1,2))./3600;
prn04_ecef = prn04_sp3(:,4:6)*1000; % [meters]
nist_ecef = [-1288398.567 -4721696.932 4078625.350]; % [meters]

[AZ_prn04, EL_prn04, RANGE_prn04] = compute_azelrange(nist_ecef, prn04_ecef);

figure(1)
plot(time, EL_prn04, 'LineWidth', 2)
grid on
xlim([0, 24])
ylim([10, 90])
ylabel("Elevation [deg]")
xlabel("Time of DOY 231 [hr]")
title("PRN04 Elevation from NIST")


%% Problem 3b

prn05_sp3 = sp3_data(sp3_data(:, 3) == 5, :);
prn05_ecef = prn05_sp3(:,4:6)*1000; % [meters]

[AZ_prn05, EL_prn05, RANGE_prn05] = compute_azelrange(nist_ecef, prn05_ecef);

figure(2)
subplot(3,1,1)
plot(time, AZ_prn05, 'LineWidth', 2)
xlim([0, 24])
ylabel("Azimuth [deg]")
grid on

subplot(3,1,2)
plot(time, EL_prn05, 'LineWidth', 2)
ylabel("Elevation [deg]")
xlim([0, 24])
grid on

subplot(3,1,3)
plot(time, RANGE_prn05, 'LineWidth', 2)
ylabel("Range [m]")
grid on
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 AER from NIST")

[~, ~, range_simple] = compute_azelrange([6378, 0, 0]*1000, [26560, 0, 0]*1000);

%% Problem 3c

prn05_visible = EL_prn05 > 15; % Decided that when elevation is greater than 15deg, sats are visible
AZ_prn05_visible = AZ_prn05(prn05_visible);
EL_prn05_visible = EL_prn05(prn05_visible);
RANGE_prn05_visible = RANGE_prn05(prn05_visible);
time_prn05_visible = time(prn05_visible);

figure(3)
subplot(3,1,1)
plot(time_prn05_visible, AZ_prn05_visible, 'o', 'LineWidth', 2)
ylabel("Azimuth [deg]")
xlim([0, 24])
grid on

subplot(3,1,2)
plot(time_prn05_visible, EL_prn05_visible, 'o', 'LineWidth', 2)
ylabel("Elevation [deg]")
xlim([0, 24])
grid on

subplot(3,1,3)
plot(time_prn05_visible, RANGE_prn05_visible, 'o', 'LineWidth', 2)
ylabel("Range [m]")
xlim([0, 24])
grid on
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 AER only when visible from NIST")

%% Problem 3d

figure(4)
subplot(3,1,1)
plot(time, AZ_prn05, 'LineWidth', 2)
hold on
plot(time, AZ_prn04, 'LineWidth', 2)
ylabel("Azimuth [deg]")
legend("PRN05", "PRN04")
grid on
xlim([0, 24])

subplot(3,1,2)
plot(time, EL_prn05, 'LineWidth', 2)
hold on
plot(time, EL_prn04, 'LineWidth', 2)
ylabel("Elevation [deg]")
grid on
xlim([0, 24])

subplot(3,1,3)
plot(time, RANGE_prn05, 'LineWidth', 2)
hold on
plot(time, RANGE_prn04, 'LineWidth', 2)
ylabel("Range [m]")
grid on
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 and PRN04 AER from NIST")

prn04_visible = EL_prn04 > 15; % Decided that when elevation is greater than 35deg, sat is visible
AZ_prn04_visible = AZ_prn04(prn04_visible);
EL_prn04_visible = EL_prn04(prn04_visible);
RANGE_prn04_visible = RANGE_prn04(prn04_visible);
time_prn04_visible = time(prn04_visible);

figure(5)
subplot(3,1,1)
plot(time_prn05_visible, AZ_prn05_visible, 'o', 'LineWidth', 2)
hold on
plot(time_prn04_visible, AZ_prn04_visible, 'o', 'LineWidth', 2)
ylabel("Azimuth [deg]")
legend("PRN05", "PRN04")
xlim([0, 24])
grid on

subplot(3,1,2)
plot(time_prn05_visible, EL_prn05_visible, 'o', 'LineWidth', 2)
hold on
plot(time_prn04_visible, EL_prn04_visible, 'o','LineWidth', 2)
ylabel("Elevation [deg]")
xlim([0, 24])
grid on

subplot(3,1,3)
plot(time_prn05_visible, RANGE_prn05_visible, 'o', 'LineWidth', 2)
hold on
plot(time_prn04_visible, RANGE_prn04_visible, 'o', 'LineWidth', 2)
ylabel("Range [m]")
xlim([0, 24])
grid on
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 and PRN04 AER only when visible from NIST")


%% Remove path

rmpath("HW2_Code/")
rmpath("HW2_DATA/")

