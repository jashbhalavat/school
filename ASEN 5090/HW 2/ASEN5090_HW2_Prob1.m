clear; clc; close all;

%% ASEN 5090 - HW 2
% Jash Bhalavat
% 09/06/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

%% Add path

addpath("HW2_Code/")
addpath("HW2_DATA/")

%% Problem 1a

filename = 'IGS0OPSFIN_20262310000_01D_15M_ORB.SP3';
sp3_data = read_sp3(filename);

prn05_sp3 = sp3_data(sp3_data(:, 3) == 5, :); % [km]
time = (prn05_sp3(:,2) - prn05_sp3(1,2))./3600; % [hr]

figure(1);
subplot(3,1,1)
plot(time, prn05_sp3(:,4)*1000, 'LineWidth', 2)
ylabel("X [m]")
xlim([0, 24])
grid on

subplot(3,1,2)
plot(time, prn05_sp3(:,5)*1000, 'LineWidth', 2)
ylabel("Y [m]")
xlim([0, 24])
grid on

subplot(3,1,3)
plot(time, prn05_sp3(:,6)*1000, 'LineWidth', 2)
ylabel("Z [m]")
xlim([0, 24])
grid on
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 ECEF Positions from SP3 file")

%% Problem 1b

filename = 'YUMA231.alm.txt';
ephem_all = read_GPSyuma(filename);
week_number = 384 * ones([length(time), 1]);
t_input =  [week_number, prn05_sp3(:,2)];
[prn05_alm_health, prn05_alm_pos, ~] = alm2pos(ephem_all, t_input, 5);

figure(2);
subplot(3,1,1)
plot(time, prn05_sp3(:,4)*1000, 'LineWidth', 2)
hold on
grid on
xlim([0, 24])
plot(time, prn05_alm_pos(:,1), '--', 'LineWidth', 2)
ylabel("X [m]")
legend("SP3", "YUMA")

subplot(3,1,2)
plot(time, prn05_sp3(:,5)*1000, 'LineWidth', 2)
hold on
grid on
xlim([0, 24])
plot(time, prn05_alm_pos(:,2), '--', 'LineWidth', 2)
ylabel("Y [m]")

subplot(3,1,3)
plot(time, prn05_sp3(:,6)*1000, 'LineWidth', 2)
hold on
grid on
xlim([0, 24])
plot(time, prn05_alm_pos(:,3), '--', 'LineWidth', 2)
ylabel("Z [m]")
xlabel("Time of DOY 231 [hr]")
sgtitle("PRN05 ECEF Positions based on SP3 and YUMA almanac")

%% Problem 1c

figure(3);
subplot(3,1,1)
plot(time, prn05_alm_pos(:,1) - prn05_sp3(:,4)*1000, 'LineWidth', 2)
ylabel("X [m]")
xlim([0, 24])
grid on

subplot(3,1,2)
plot(time, prn05_alm_pos(:,2) - prn05_sp3(:,5)*1000, 'LineWidth', 2)
ylabel("Y [m]")
xlim([0, 24])
grid on

subplot(3,1,3)
plot(time, prn05_alm_pos(:,3) - prn05_sp3(:,6)*1000, 'LineWidth', 2)
ylabel("Z [m]")
xlabel("Time of DOY 231 [hr]")
grid on
xlim([0, 24])
sgtitle("PRN05 ECEF Position difference between almanac and SP3")


%% Remove path

rmpath("HW2_Code/")
rmpath("HW2_DATA/")