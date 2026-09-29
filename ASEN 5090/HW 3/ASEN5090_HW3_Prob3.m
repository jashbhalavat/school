clear; clc; close all;

%% ASEN 5090 - HW 3
% Jash Bhalavat
% 09/20/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1
% Earth rotation rate
omega_dot = 7.2921151467e-5; % rad/s, Lecture 5, slide 31

nist_ecef = [-1288398.567 -4721696.932 4078625.350]; % [meters]

%% Add path

addpath("HW3_Code/")
addpath("HW3_DATA/")
addpath("HW2_Code/")
addpath("HW2_DATA/")
addpath("../HW 2")

%% Problem 3a

% Read broadcast file and save data into array
broadcast_ephem_filename = 'brdc2310.26n';
ephem_data = read_clean_GPSbroadcast(broadcast_ephem_filename);

% Set filename and use rinexread to extract GPS data
nist_rinex_filename = 'NIST00USA_R_20262310000_01D_30S_MO.rnx';
rinex_data = rinexread(nist_rinex_filename).GPS;

% Get time and convert data to arrays
time = rinex_data.Time;
rinex_data_array = table2array(rinex_data);

% Only extract PRN05 data and time
prn05_obs_data = rinex_data_array(rinex_data_array(:,1) == 5, :);
prn05_obs_time = time(rinex_data_array(:,1) == 5);
prn05_obs_time.TimeZone = 'UTC';

% Define the absolute GPS Epoch (January 6, 1980 is a Sunday)
gpsEpoch = datetime(1980, 1, 6, 'TimeZone', 'UTC');
gpsEpochArray = repmat(gpsEpoch, size(prn05_obs_time));

% Calculate total days elapsed since the GPS epoch
totalDays = floor(days(prn05_obs_time - gpsEpochArray));

% Find the GPS Week Number 
gpsWeek = floor(totalDays / 7);

% Find the start of the current GPS week (the previous Sunday at 00:00:00)
currentSunday = gpsEpoch + days(gpsWeek * 7);

% Calculate Time of Week (TOW) in seconds
TOW = seconds(prn05_obs_time - currentSunday);

week_number = 2432 * ones([length(TOW), 1]);
t_input =  [week_number, TOW];

[~, prn05_broadcast] = eph2pvt2025(ephem_data, t_input, 5);

[prn05_az_R0, prn05_el_R0, prn05_range_R0] = compute_azelrange(nist_ecef, prn05_broadcast);

plot_time = (t_input(:,2)-t_input(1,2))/3600;

figure(1)
subplot(3,1,1)
plot(plot_time, prn05_az_R0,'o','MarkerSize',3)
grid on
ylabel("Azimuth [deg]")
xlim([0, 24])

subplot(3,1,2)
plot(plot_time, prn05_el_R0,'o','MarkerSize',3)
grid on
xlim([0, 24])
ylabel("Elevation [deg]")

subplot(3,1,3)
plot(plot_time, prn05_range_R0,'o','MarkerSize',3)
grid on
ylabel("Range [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 Az-El-Range (R0)")

%% Problem 3b

prn05_range_R1 = compute_expected_range(t_input, prn05_broadcast, nist_ecef, c, ephem_data, 5, omega_dot);

figure(2)
subplot(1,2,1)
plot(plot_time, prn05_range_R0,'o','MarkerSize',3)
hold on
plot(plot_time, prn05_range_R1,'o','MarkerSize',3)
grid on
ylabel("Range [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
legend("R0", "R1")
title("R0 and R1 ranges")

subplot(1,2,2)
plot(plot_time, prn05_range_R1 - prn05_range_R0,'o','MarkerSize',3)
grid on
ylabel("Range difference [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
title("R1 minus R0")

largest_difference = max(abs(prn05_range_R1 - prn05_range_R0));

%% Problem 4a

figure(3)
subplot(1,2,1)
plot(plot_time, prn05_obs_data(:,6), 'o', 'MarkerSize', 3)
hold on
plot(plot_time, prn05_range_R1, 'o', 'MarkerSize', 3)
grid on
xlim([0, 24])
ylabel("Ranges [m]")
xlabel("Time of DOY 231 [hr]")
legend("C1C", "R1")
title("Observed and Computed Range")
sgtitle("PRN05 Observed vs Computed Range")

subplot(1,2,2)
plot(plot_time, prn05_obs_data(:,6) - prn05_range_R1, 'o', 'LineWidth',0.25)
xlabel("Time")
xlim([0, 24])
ylabel("C1C Pseudorange [meters]")
grid on
title("Observed minus Computed vs Time")


%% Remove path

rmpath("HW2_Code/")
rmpath("HW2_DATA/")
rmpath("HW3_Code/")
rmpath("HW3_DATA/")
rmpath("../HW 2")
