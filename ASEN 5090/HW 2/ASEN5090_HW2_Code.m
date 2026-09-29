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

clear; clc; close all;

%% ASEN 5090 - HW 2
% Jash Bhalavat
% 09/10/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

%% Problem 2a

nist_ecef = [-1288398.567 -4721696.932 4078625.350]; % [meters]
smead_lla = [40.010331, -105.244285, 1600]; % [deg, deg, meters]

nist_lla = ecef2lla(nist_ecef); % [deg, deg, meters]
smead_ecef = lla2ecef(smead_lla); % [meters]

equa_lla = [0, nist_lla(2), 10]; % [deg, deg, meters]
equa_ecef = lla2ecef(equa_lla); % [meters]

%% Problem 2b

r_hat_ecef_nist = nist_ecef./norm(nist_ecef);
C_ECEF2ENU_NIST = ECEF2ENU(nist_lla(1), nist_lla(2));
r_hat_enu_nist = C_ECEF2ENU_NIST * r_hat_ecef_nist';

r_hat_ecef_equa = equa_ecef./norm(equa_ecef);
C_ECEF2ENU_EQUA = ECEF2ENU(equa_lla(1), equa_lla(2));
r_hat_enu_equa = C_ECEF2ENU_EQUA * r_hat_ecef_equa';

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

%% User ECEF (problem 2a)

nist_ecef = [-1288398.567 -4721696.932 4078625.350]; % [meters]
smead_lla = [40.010331, -105.244285, 1600]; % [deg, deg, meters]

nist_lla = ecef2lla(nist_ecef); % [deg, deg, meters]
smead_ecef = lla2ecef(smead_lla); % [meters]

equa_lla = [0, nist_lla(2), 10]; % [deg, deg, meters]
equa_ecef = lla2ecef(equa_lla); % [meters]

%% Problem 4a

filename = 'IGS0OPSFIN_20262310000_01D_15M_ORB.SP3';
sp3_data = read_sp3(filename);

all_sats_ecef = sp3_data(:,4:6)*1000; % [meters]

[AZ_NIST, EL_NIST, RANGE_NIST] = compute_azelrange(nist_ecef, all_sats_ecef);

all_sats_visible_NIST = EL_NIST > 35;
AZ_NIST_visible = AZ_NIST(all_sats_visible_NIST);
EL_NIST_visible = EL_NIST(all_sats_visible_NIST);
prn_NIST_visible = sp3_data(all_sats_visible_NIST, 3);

figure(1)
plotAzEl(AZ_NIST_visible, EL_NIST_visible, prn_NIST_visible)
title("NIST Skyplot All Sats")

%% Problem 4b

% SMEAD visibility
[AZ_SMEAD, EL_SMEAD, RANGE_SMEAD] = compute_azelrange(smead_ecef, all_sats_ecef);

all_sats_visible_SMEAD = EL_SMEAD > 10;
AZ_SMEAD_visible = AZ_SMEAD(all_sats_visible_SMEAD);
EL_SMEAD_visible = EL_SMEAD(all_sats_visible_SMEAD);
prn_SMEAD_visible = sp3_data(all_sats_visible_SMEAD, 3);

figure(2)
plotAzEl(AZ_SMEAD_visible, EL_SMEAD_visible, prn_SMEAD_visible)
title("SMEAD Skyplot All Sats")

% EQUA visibility
[AZ_EQUA, EL_EQUA, RANGE_EQUA] = compute_azelrange(equa_ecef, all_sats_ecef);

all_sats_visible_EQUA = EL_EQUA > 10;
AZ_EQUA_visible = AZ_EQUA(all_sats_visible_EQUA);
EL_EQUA_visible = EL_EQUA(all_sats_visible_EQUA);
prn_EQUA_visible = sp3_data(all_sats_visible_EQUA, 3);

figure(3)
plotAzEl(AZ_EQUA_visible, EL_EQUA_visible, prn_EQUA_visible)
title("EQUA Skyplot All Sats")


%% Remove path

rmpath("HW2_Code/")
rmpath("HW2_DATA/")

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

