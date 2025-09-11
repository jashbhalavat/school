clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 3, Problem 1
% Jash Bhalavat

%% Given

r_c_N = [-4893.268,3864.478,3169.646]'; % km
v_c_N = [-3.91386, -6.257673, 1.59797]'; % km/s

r_d_N = [-4892.98,3863.073,3170.619]'; % km
v_d_N = [-3.913302,-6.258661,1.598199]'; % km/s

mu = 398600; % km3/s2

[rho_O, rho_O_prime] = inertial_cd_2_hill_rel(r_c_N, v_c_N, r_d_N, v_d_N);

% # adjust the return matrix values as needed
% def result():
%     rho_H = [-0.536809605885212, 1.22119508201654, 1.10644295772025]  # km
%     rhoP_H = [0.000486025146410436, 0.00115805696520614, 0.000558863421179336] # km/s
%     return rho_H, rhoP_H

function [rho_O, rho_O_prime] = inertial_cd_2_hill_rel(r_c_N, v_c_N, r_d_N, v_d_N)
   ohat_r = r_c_N/norm(r_c_N);
   h_N = cross(r_c_N, v_c_N);
   ohat_h = h_N/norm(h_N);
   ohat_theta = cross(ohat_h, ohat_r);

   dcm_ON = [ohat_r'; ohat_theta'; ohat_h'];

   fdot = norm(h_N)/(norm(r_c_N))^2;
   omega_ON_O = [0; 0; fdot];

   rho_O = dcm_ON * (r_d_N - r_c_N);
   rho_O_prime = dcm_ON * (v_d_N - v_c_N) - cross(omega_ON_O, rho_O);
end
