clear; clc; close all;

%% ASEN 5090 - HW 3
% Jash Bhalavat
% 09/20/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

%% Add path

addpath("HW3_Code/")
addpath("HW3_DATA/")
addpath("HW2_Code/")
addpath("HW2_DATA/")

%% Problem 2a

% Read broadcast file and save data into array
broadcast_ephem_filename = 'brdc2310.26n';
ephem_data = read_clean_GPSbroadcast(broadcast_ephem_filename);

%% Problem 2b

filename = 'IGS0OPSFIN_20262310000_01D_15M_ORB.SP3';
sp3_data = read_sp3(filename);

prn05_sp3 = sp3_data(sp3_data(:, 3) == 5, :); % [km]
time = (prn05_sp3(:,2) - prn05_sp3(1,2))./3600; % [hr]

week_number = 2432 * ones([length(time), 1]);
t_input =  [week_number, prn05_sp3(:,2)];

[~, prn05_broadcast, ~, prn05_clock_bias] = eph2pvt2025(ephem_data, t_input, 5);

figure(1)
subplot(3,2,1)
plot(time, prn05_sp3(:,4)*1000, 'LineWidth',2)
hold on
plot(time, prn05_broadcast(:,1), 'LineWidth',2)
grid on
ylabel("X [m]")
xlim([0, 24])
legend("SP3", "Broadcast Ephemeris")
title("SP3 and Broadcast Ephemeris Positions")

subplot(3,2,3)
plot(time, prn05_sp3(:,5)*1000, 'LineWidth',2)
hold on
plot(time, prn05_broadcast(:,2), 'LineWidth',2)
grid on
xlim([0, 24])
ylabel("Y [m]")

subplot(3,2,5)
plot(time, prn05_sp3(:,6)*1000, 'LineWidth',2)
hold on
plot(time, prn05_broadcast(:,3), 'LineWidth',2)
grid on
ylabel("Z [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 Position from Broadcast Ephmeris vs SP3")

subplot(3,2,2)
plot(time, prn05_sp3(:,4)*1000 - prn05_broadcast(:,1), 'LineWidth',2)
grid on
ylabel("X [m]")
xlim([0, 24])
title("SP3 minus Broadcast Ephemeris Positions")

subplot(3,2,4)
plot(time, prn05_sp3(:,5)*1000 - prn05_broadcast(:,2), 'LineWidth',2)
grid on
ylabel("Y [m]")
xlim([0, 24])

subplot(3,2,6)
plot(time, prn05_sp3(:,6)*1000 - prn05_broadcast(:,3), 'LineWidth',2)
grid on
ylabel("Z [m]")
xlabel("Time of DOY 231 [hr]")

%% Problem 2c

figure(2)
plot(time, prn05_clock_bias, 'LineWidth',2)
grid on
ylabel("Satellite Clock Bias [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
title("PRN05 Clock Bias from Broadcast Ephemeris")


%% Remove path

rmpath("HW2_Code/")
rmpath("HW2_DATA/")
rmpath("HW3_Code/")
rmpath("HW3_DATA/")
