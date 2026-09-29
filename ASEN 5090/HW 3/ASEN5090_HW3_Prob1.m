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

%% Problem 1a

% Set filename and use rinexread to extract GPS data
nist_rinex_filename = 'NIST00USA_R_20262310000_01D_30S_MO.rnx';
rinex_data = rinexread(nist_rinex_filename).GPS;

% Get time and convert data to arrays
time = rinex_data.Time;
rinex_data_array = table2array(rinex_data);

% Only extract PRN05 data and time
prn05_data = rinex_data_array(rinex_data_array(:,1) == 5, :);
prn05_time = time(rinex_data_array(:,1) == 5);

%% Problem 1b

figure(1)
plot(prn05_time, prn05_data(:,6), 'o', 'LineWidth',0.25)
xlabel("Time")
ylabel("C1C Pseudorange [meters]")
grid on
title("PRN05 C1C pseudorange vs Time")

%% Problem 1c

figure(2)
plot(prn05_time, prn05_data(:,13), 'o', 'LineWidth',0.25)
xlabel("Time")
ylabel("S1C Signal to Noise Ratio [dB-Hz]")
grid on
title("PRN05 S1C Signal to Noise Ratio vs Time")

%% Problem 1d

prn05_l1c_cycles = prn05_data(:,8); % [cycles]
l1_frequency = 1575.42e6; % [hz]
l1_wavelength = c / l1_frequency; % [meters]
prn05_l1c_meters = prn05_l1c_cycles * l1_wavelength;

figure(3)
subplot(2, 2, [1,3])
plot(prn05_time, prn05_data(:,6), 'o', 'LineWidth',0.25)
hold on
plot(prn05_time, prn05_data(:,28), 'o', 'LineWidth',0.25)
plot(prn05_time, prn05_l1c_meters, 'o', 'LineWidth',0.25)
legend("C1C", "C2L", "L1C")
xlabel("Time")
ylabel("Pseudorange [m]")
grid on
title("C1C, C2L, L1C (converted to meters) vs Time")

subplot(2,2,2)
plot(prn05_time, prn05_data(:,6)-prn05_data(:,28))
grid on
title("C1C minus C2L vs Time")
ylabel("Pseudorange [m]")
xlabel("Time")

subplot(2,2,4)
plot(prn05_time, prn05_data(:,6)-prn05_l1c_meters)
grid on
title("C1C minus L1C (converted to meters) vs Time")
ylabel("Pseudorange [m]")
xlabel("Time")
sgtitle("PRN05")

%% Remove path

rmpath("HW3_Code/")
rmpath("HW3_DATA/")
