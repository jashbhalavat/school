clear; clc; close all;

%% ASEN 5090 - HW 4
% Jash Bhalavat
% 09/25/2026

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

%% Save Rinex Data

% Set filename and use rinexread to extract GPS data
nist_rinex_filename = 'NIST00USA_R_20262310000_01D_30S_MO.rnx';
rinex_data = rinexread(nist_rinex_filename).GPS;

save("rinex_data.mat", 'rinex_data')


%% Problem 1

load('rinex_data.mat')

% Get time and convert data to arrays
time = rinex_data.Time;
rinex_data_array = table2array(rinex_data);

% Only extract PRN05 data and time
prn14_data = rinex_data_array(rinex_data_array(:,1) == 14, :);
% Kill first observation
prn14_data = prn14_data(2:end,:);
prn14_C1C = prn14_data(:,6); % meters
prn14_time = time(rinex_data_array(:,1) == 14);
prn14_time = prn14_time(2:end,:);
prn14_time.TimeZone = 'UTC';

% Read broadcast file and save data into array
broadcast_ephem_filename = 'brdc2310.26n';
ephem_data = read_clean_GPSbroadcast(broadcast_ephem_filename);

% Define the absolute GPS Epoch (January 6, 1980 is a Sunday)
gpsEpoch = datetime(1980, 1, 6, 'TimeZone', 'UTC');
gpsEpochArray = repmat(gpsEpoch, size(prn14_time));

% Calculate total days elapsed since the GPS epoch
totalDays = floor(days(prn14_time - gpsEpochArray));

% Find the GPS Week Number 
gpsWeek = floor(totalDays / 7);

% Find the start of the current GPS week (the previous Sunday at 00:00:00)
currentSunday = gpsEpoch + days(gpsWeek * 7);

% Calculate Time of Week (TOW) in seconds
TOW = seconds(prn14_time - currentSunday);

week_number = 2432 * ones([length(TOW), 1]);
t_input =  [week_number, TOW];

[~, prn14_range_R0] = eph2pvt2025(ephem_data, t_input, 14);

plot_time = (t_input(:,2)-t_input(1,2))/3600;

prn14_range_R1 = compute_expected_range(t_input, prn14_range_R0, nist_ecef, c, ephem_data, 14, omega_dot);

dPR0 = prn14_C1C - prn14_range_R1;

figure(1)
plot(prn14_time, dPR0,'o','MarkerSize',3)
grid on
ylabel("Range [m]")
xlabel("Time of DOY 231")
title("dPR0")

%% Problem 2

[~, prn14_range_R0, ~, prn14_bsv] = eph2pvt2025(ephem_data, t_input, 14);

dPR1 = prn14_C1C - (prn14_range_R1 - prn14_bsv);

figure(2)
plot(prn14_time, dPR1,'o','MarkerSize',3)
grid on
ylabel("Range [m]")
xlabel("Time of DOY 231")
title("dPR1")

%% Problem 3

[~, prn14_range_R0, ~, prn14_bsv, ~, ~, prn14_relsv] = eph2pvt2025(ephem_data, t_input, 14);

dPR2 = prn14_C1C - (prn14_range_R1 - prn14_bsv - prn14_relsv);

figure(2)
plot(prn14_time, dPR2,'o','MarkerSize',3)
grid on
ylabel("Range [m]")
xlabel("Time of DOY 231")
title("dPR2")

%% Problem 4




%% Functions

function [tropo] = tropomodel(zd, elevation)


end


%% Remove path

rmpath("HW2_Code/")
rmpath("HW2_DATA/")
rmpath("HW3_Code/")
rmpath("HW3_DATA/")
rmpath("../HW 2")