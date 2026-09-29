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
