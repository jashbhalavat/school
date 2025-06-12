clear; clc; close all;

% ASEN 5010 - Project, Programming Assignment, Task 2
% Spring 2025
% Jash Bhalavat

%% Given

r_mars = 3396.19; % km
h = 400; % km
r_lmo = r_mars + h;

mu_mars = 42828.3; % km^3/s^2
theta_dot_lmo = sqrt(mu_mars/r_lmo^3);

r_gmo = 20424.2; % km
theta_dot_gmo = 0.0000709003; % rad/s

antenna_n_hat_b = [-1, 0, 0]';
sensor_n_hat_b = [1, 0, 0]';
sa_n_hat_b = [0, 0, 1]';

sigma_BN_t0 = [0.3, -0.4, 0.5]';
omega_BN_B_t0 = [1, 1.75, -2.2]' * pi/180; % rad/s
I_B = diag([10, 5, 7.5]); % kg-m^2

omega_0_lmo = 20 * pi/180;
i_0_lmo = 30 * pi/180;
theta_0_lmo = 60 * pi/180;
period_lmo = 2*pi*sqrt(r_lmo^3/mu_mars);

omega_0_gmo = 0 * pi/180;
i_0_gmo = 0 * pi/180;
theta_0_gmo = 250 * pi/180;
period_gmo = 2*pi*sqrt(r_gmo^3/mu_mars);

n_lmo = ceil(period_lmo); % time steps per orbit
n_gmo = ceil(period_gmo);

t_lmo = 0:n_lmo;
t_gmo = 0:n_gmo;

%% Task 1

for i = 1:length(t_lmo)
    theta_lmo(i,1) = wrapTo2Pi(theta_dot_lmo*t_lmo(i) + theta_0_lmo);
    ea_ON_lmo(i,:) = [omega_0_lmo, i_0_lmo, theta_lmo(i,1)];
    [r_N_lmo(i,:), r_dot_N_lmo(i,:)] = inertial_velocity(r_lmo, ea_ON_lmo(i,:), theta_dot_lmo);
end

for i = 1:length(t_gmo)
    theta_gmo(i,1) = wrapTo2Pi(theta_dot_gmo*t_gmo(i) + theta_0_gmo);
    ea_ON_gmo(i,:) = [omega_0_gmo, i_0_gmo, theta_gmo(i,1)];
    [r_N_gmo(i,:), r_dot_N_gmo(i,:)] = inertial_velocity(r_gmo, ea_ON_gmo(i,:), theta_dot_gmo);
end


%% Task 2

% [NT] = [t1_N, t2_N, t3_N]
% Therefore, [NH] = [ir_N, itheta_N, ih_N];

% @300 seconds, t_lmo index = 301
t_300_ind = 301;

DCM_HN_at_300 = orbit_frame_orientation(t_300_ind, r_N_lmo, r_dot_N_lmo);

%% Functions

function DCM_HN = orbit_frame_orientation(t, r_vec, r_dot_vec)
    r = r_vec(t,:)';
    r_dot = r_dot_vec(t,:)';
    i_r_N = r / norm(r);
    i_h_N = cross(r, r_dot) / norm(cross(r, r_dot));
    i_theta_N = cross(i_h_N, i_r_N);
    NH = [i_r_N, i_theta_N, i_h_N];
    DCM_HN = NH';
end


function [r_N, r_dot_N] = inertial_velocity(r, ea, theta_dot)
    % Inertial s/c velocity vector for a circular orbit
    % Inputs
    % r - radius
    % ea - 3-1-1 euler angles
    % 
    % Outputs
    % r_N - inertial position
    % r_dot_N - inertial velocity

    r_O = [r, 0, 0]';
    omega_ON_O = [0, 0, theta_dot]';
    
    ON = R3(ea(3)) * R1(ea(2)) * R3(ea(1));


    r_N = ON' * r_O;
    r_dot_N = ON' * cross(omega_ON_O, r_O);

end

function state_dot = diff_eq(state, I)
    % 6x1 state EOMs
    % Inputs:
    % state - [sigma_BN'; omega_BN_B']
    % I - 3x3 principal MOI [kg-m^2]
    % 
    % Output:
    % state_dot - [sigma_BN_dot', omega_BN_B_dot']
    sigma_BN = state(1:3);
    omega_BN_B = state(4:6);

    sigma_sq = norm(sigma_BN)^2;

    omega_BN_B_dot(1,1) = -(I(3,3) - I(2,2))/I(1,1) * omega_BN_B(2) * omega_BN_B(3);
    omega_BN_B_dot(2,1) = -(I(1,1) - I(3,3))/I(2,2) * omega_BN_B(1) * omega_BN_B(3);
    omega_BN_B_dot(3,1) = -(I(2,2) - I(1,1))/I(3,3) * omega_BN_B(2) * omega_BN_B(1);

    sigma_BN_dot = 1/4 * ((1-sigma_sq)*eye(3) + 2*skew_symmetric(sigma) + 2*sigma*sigma') * omega_BN_B;

    state_dot = [sigma_BN_dot; omega_BN_B_dot];
end


