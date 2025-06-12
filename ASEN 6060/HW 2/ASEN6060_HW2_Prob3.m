clear; clc; close all;

% ASEN 6060 - HW 2, Prob 3
% Jash Bhalavat
% Spring 2025

mu = linspace(1e-7, 0.5, 1000);

for i = 1:length(mu)

    % Find L4 and L5 points
    x_L4 = 1/2 - mu(i);
    y_L4 = sqrt(3)/2;
    y_l5 = -sqrt(3)/2;
    x_Ls = [x_L4, y_L4; x_L4, y_l5];
    
    % Calculate out of plane modes for all L4
    eq_pts_in_plane_modes_l4(i,:) = in_plane_modes(mu(i), x_Ls(1,:));

    % Calculate out of plane modes for all L5
    eq_pts_in_plane_modes_l5(i,:) = in_plane_modes(mu(i), x_Ls(2,:));


end
