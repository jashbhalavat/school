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