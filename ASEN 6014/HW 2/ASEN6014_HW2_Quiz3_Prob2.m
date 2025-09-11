clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 3, Problem 2
% Jash Bhalavat

%% Given

r_c_N = [-4893.268,3864.478,3169.646]'; % km
v_c_N = [-3.91386, -6.257673, 1.59797]'; % km/s

rho_O = [-0.537,1.221,1.106]'; % km
rho_O_prime = [0.000486,0.001158,0.0005590]'; % km/s

mu = 398600; % km3/s2

[r_d_N, v_d_N] = hill_rel_2_inertial_deputy(r_c_N, v_c_N, rho_O, rho_O_prime);

% # adjust the return matrix values as needed
% def result():
%     rd_N = [-4892.97998390252, 3863.07309498812, 3170.61848888692]  # km
%     vd_N = [-3.91330192656190, -6.25866069175413, 1.59819914685574] # km/s
%     return rd_N, vd_N



function [r_d_N, v_d_N] = hill_rel_2_inertial_deputy(r_c_N, v_c_N, rho_O, rho_O_prime)
   ohat_r = r_c_N/norm(r_c_N);
   h_N = cross(r_c_N, v_c_N);
   ohat_h = h_N/norm(h_N);
   ohat_theta = cross(ohat_h, ohat_r);

   dcm_ON = [ohat_r'; ohat_theta'; ohat_h'];
   dcm_NO = dcm_ON';

   fdot = norm(h_N)/(norm(r_c_N))^2;
   omega_ON_O = [0; 0; fdot];

   r_d_N = r_c_N + dcm_NO*rho_O;
   v_d_N = v_c_N + dcm_NO * (rho_O_prime + cross(omega_ON_O, rho_O));
end
