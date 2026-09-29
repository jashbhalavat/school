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

%% Problem 5

% Read broadcast file and save data into array
broadcast_ephem_filename = 'brdc2310.26n';
ephem_data = read_clean_GPSbroadcast(broadcast_ephem_filename);

% Set filename and use rinexread to extract GPS data
nist_rinex_filename = 'NIST00USA_R_20262310000_01D_30S_MO.rnx';
rinex_data = rinexread(nist_rinex_filename).GPS;

% Get time and convert data to arrays
time = rinex_data.Time;
rinex_data_array = table2array(rinex_data);

% Only extract PRN20 data and time
prn20_obs_data = rinex_data_array(rinex_data_array(:,1) == 20, :);
prn20_obs_time = time(rinex_data_array(:,1) == 20);
prn20_obs_time.TimeZone = 'UTC';

% Define the absolute GPS Epoch (January 6, 1980 is a Sunday)
gpsEpoch = datetime(1980, 1, 6, 'TimeZone', 'UTC');
gpsEpochArray = repmat(gpsEpoch, size(prn20_obs_time));

% Calculate total days elapsed since the GPS epoch
totalDays = floor(days(prn20_obs_time - gpsEpochArray));

% Find the GPS Week Number 
gpsWeek = floor(totalDays / 7);

% Find the start of the current GPS week (the previous Sunday at 00:00:00)
currentSunday = gpsEpoch + days(gpsWeek * 7);

% Calculate Time of Week (TOW) in seconds
TOW = seconds(prn20_obs_time - currentSunday);

week_number = 2432 * ones([length(TOW), 1]);
t_input =  [week_number, TOW];

[~, prn20_broadcast] = eph2pvt2025(ephem_data, t_input, 20);

[prn20_az_R0, prn20_el_R0, prn20_range_R0] = compute_azelrange(nist_ecef, prn20_broadcast);

plot_time = (t_input(:,2)-t_input(1,2))/3600;

figure(1)
subplot(3,1,1)
plot(plot_time, prn20_az_R0,'o','MarkerSize',3)
grid on
ylabel("Azimuth [deg]")
xlim([0, 24])

subplot(3,1,2)
plot(plot_time, prn20_el_R0,'o','MarkerSize',3)
grid on
xlim([0, 24])
ylabel("Elevation [deg]")

subplot(3,1,3)
plot(plot_time, prn20_range_R0,'o','MarkerSize',3)
grid on
ylabel("Range [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN20 Az-El-Range (R0)")

prn20_range_R1 = compute_expected_range(t_input, prn20_broadcast, nist_ecef, c, ephem_data, 20, omega_dot);

figure(2)
subplot(1,2,1)
plot(plot_time, prn20_range_R0,'o','MarkerSize',3)
hold on
plot(plot_time, prn20_range_R1,'o','MarkerSize',3)
grid on
ylabel("Range [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
legend("R0", "R1")
title("R0 and R1 ranges")

subplot(1,2,2)
plot(plot_time, prn20_range_R1 - prn20_range_R0,'o','MarkerSize',3)
grid on
ylabel("Range difference [m]")
xlim([0, 24])
xlabel("Time of DOY 231 [hr]")
title("R1 minus R0")

largest_difference = max(abs(prn20_range_R1 - prn20_range_R0));

figure(3)
subplot(1,2,1)
plot(plot_time, prn20_obs_data(:,6), 'o', 'MarkerSize', 3)
hold on
plot(plot_time, prn20_range_R1, 'o', 'MarkerSize', 3)
grid on
xlim([0, 24])
ylabel("Ranges [m]")
xlabel("Time of DOY 231 [hr]")
legend("C1C", "R1")
title("Observed and Computed Range")
sgtitle("PRN20 Observed vs Computed Range")

subplot(1,2,2)
plot(plot_time, prn20_obs_data(:,6) - prn20_range_R1, 'o', 'LineWidth',0.25)
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


%% FUNCTIONS

function R1 = compute_expected_range(t_r, r_gps_tr, r_rx, c, ephem_data, prn, omega_E)
    % Computed expected range given Tr [WN, TOW], r_gps [ECEF, m], c [m/s]
    max_count = 5;

    % Initial measure of time taken for signal to get from satellite to
    % receiver
    dt = norm(r_gps_tr(1,:) - r_rx)/c;

    for j = 1:max_count
        for i = 1:length(t_r)
            % Step 2
            R(i) = norm(r_gps_tr(i,:) - r_rx);
    
            % Step 3
            t_t(i,1) = t_r(i,2) - R(i)/c;
        end
        t_t = [t_r(:,1), t_t];
    
        % Step 4
        [~, r_gps_tt] = eph2pvt2025(ephem_data, t_t, prn);
    
        % Step 5
        for i = 1:length(t_r)
            phi = omega_E * (t_r(i,2)-t_t(i,2));
            r_gps_tr(i,:) = ([cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1] * r_gps_tt(i,:)')';
    
            % Step 6
            R1(i,1) = norm(r_gps_tr(i,:) - r_rx);
        end

        % Calculate new signal travel time
        dt_new = R1(1)/c;

        % Compare signal travel times and break if within tolerance
        if abs(dt_new - dt) < 1e-12
            break
        end

        % If not within tolerance, update current dt
        dt = dt_new;
    end
end

