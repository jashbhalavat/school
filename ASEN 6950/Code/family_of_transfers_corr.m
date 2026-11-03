clear; clc; close all

%% Constants and IC
    
G = 6.67408 * 10^-11; % m3/(kgs2)
G = G / (10^9); % km3/(kgs2)

% Earth
mu_earth = 398600.435507; % km3/s2
a_earth = 149598023; % km
e_earth = 0.016708617;
mass_earth = mu_earth / G; % kg

% Moon
mu_moon = 4902.800118; % km3/s2
a_moon = 384400; % km
e_moon = 0.05490;
mass_moon = mu_moon / G; % kg

% Earth-Moon system
mass_ratio_em = mass_moon / (mass_earth + mass_moon);
m_star_em = mass_earth + mass_moon; % kg
l_star_em = a_moon; % km 
t_star_em = sqrt(l_star_em^3/(G * m_star_em)); % s
mu = mass_ratio_em;

T = 15e-3/1000; % kN
Isp = 1320; % sec

init_mass = 100;

% Non-dim thrust
f = (T*t_star_em^2)/(l_star_em*init_mass);

% g is assumed to be 9.80665e-3 km/s^2
mdot = -(f*l_star_em)/(Isp*9.80665e-3*t_star_em);

l1_lyapunov_orbits = load("V_family_L1_Lyapunov.mat").V_family;
l2_lyapunov_orbits = load("V_family_L2_Lyapunov_orbits.mat").V_family;
init_orbit_idx = 25;
final_orbit_idx = 40;
em_eq_pts = load("eq_points.mat").em_eq_pts;

l2_pos = [em_eq_pts(2,:), 0];
l1_pos = [em_eq_pts(1,:), 0];
p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

function [value, isterminal, direction] = single_cross_direction_agnostic(t, y, mu)
    value = y(1) - (1-mu);    % Detect when y(1) crosses 1-mu (moon)
    direction = 0;   % Trigger for both increasing and decreasing
    isterminal = 1;  % Stop after one event
end

function [value, isterminal, direction] = multiple_cross_direction_agnostic(t, y, mu)
    value = y(1) - (1-mu);    % Detect when y(1) crosses 1-mu (moon)
    direction = 0;   % Trigger for both increasing and decreasing
    isterminal = 0;  % Don't stop at event trigger
end
    
% Set options for ode113
options_no_events = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
options_single_cross = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y)single_cross_direction_agnostic(t,y,mu));
options_mult_cross = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y)multiple_cross_direction_agnostic(t,y,mu));

[tout_lyapunov_init, xout_lyapunov_init] = ode113(@(t,state)CR3BP(state, mu), [0, l1_lyapunov_orbits(7,init_orbit_idx)], l1_lyapunov_orbits(1:6,init_orbit_idx), options_no_events);

[tout_lyapunov_final, xout_lyapunov_final] = ode113(@(t,state)CR3BP(state, mu), [0, l2_lyapunov_orbits(7,final_orbit_idx)], l2_lyapunov_orbits(1:6,final_orbit_idx), options_no_events);

lyapunov_init_jacobi = jacobiConstantCR3BP(xout_lyapunov_init, mu);
lyapunov_final_jacobi = jacobiConstantCR3BP(xout_lyapunov_final, mu);

num_angles = 50;
angles = linspace(0, 2*pi, num_angles);

% Straight down vector
neg_y_vec = [0, -1, 0]';

for i = 1:num_angles
    dcm = R3(angles(i));
    thrust_direction_50(:,i) = dcm * neg_y_vec;
end


num_angles = 25;
angles = linspace(0, 2*pi, num_angles);

for i = 1:num_angles
    dcm = R3(angles(i));
    thrust_direction_25(:,i) = dcm * neg_y_vec;
end

num_pre_moon = 20;
num_post_moon = 20;

%% Find files

load('saved_final_state_init_mult_cross.mat');
load('saved_final_state_init_single_cross.mat');

% Find all files starting with "PrefixName" in the current folder
file_list = dir('final_nondim_mass_*');

% Database of final states
final_states_single = [];
final_states_mult = [];

for i = 31:length(file_list)
    file_name = file_list(i).name;
    final_mass(i,1) = str2double(file_name(end-15:end-4));

    single_cross = load(file_name).saved_final_state_final_single_cross;
    % single_cross(:,7) = final_mass(i);
    % final_states_single = [final_states_single; single_cross];

    mult_cross = load(file_name).saved_final_state_final_mult_cross;
    % mult_cross(:,8) = final_mass(i);
    % final_states_mult = [final_states_mult; mult_cross];

    single_single = compare_poincare_maps_fnorm(saved_final_state_init_single_cross, single_cross);
    single_mult = compare_poincare_maps_fnorm(saved_final_state_init_single_cross, mult_cross);
    mult_single = compare_poincare_maps_fnorm(saved_final_state_init_mult_cross, single_cross);
    mult_mult = compare_poincare_maps_fnorm(saved_final_state_init_mult_cross, mult_cross);

    option_array = [single_single(1), single_mult(1), mult_single(1), mult_mult(1)];
    if find(option_array==min(option_array)) == 1
        uncorrected_init_idx = single_single(2);
        uncorrected_final_idx = single_single(3);
        uncorrected_init_i = saved_final_state_init_single_cross(uncorrected_init_idx,4);
        uncorrected_init_j = saved_final_state_init_single_cross(uncorrected_init_idx,5);
        uncorrected_final_i = single_cross(uncorrected_final_idx,4);
        uncorrected_final_j = single_cross(uncorrected_final_idx,5);
        % uncorrected_init_time = saved_final_state_init_single_cross(uncorrected_init_idx,6);
        % uncorrected_final_time = saved_final_state_final_single_cross(uncorrected_final_idx,6);
    elseif find(option_array==min(option_array)) == 2
        uncorrected_init_idx = single_mult(2);
        uncorrected_final_idx = single_mult(3);
        uncorrected_init_i = saved_final_state_init_single_cross(uncorrected_init_idx,4);
        uncorrected_init_j = saved_final_state_init_single_cross(uncorrected_init_idx,5);
        uncorrected_final_i = mult_cross(uncorrected_final_idx,4);
        uncorrected_final_j = mult_cross(uncorrected_final_idx,5);
        % uncorrected_init_time = saved_final_state_init_single_cross(uncorrected_init_idx,6);
        uncorrected_final_time = mult_cross(uncorrected_final_idx,6);
    elseif find(option_array==min(option_array)) == 3
        uncorrected_init_idx = mult_single(2);
        uncorrected_final_idx = mult_single(3);
        uncorrected_init_i = saved_final_state_init_mult_cross(uncorrected_init_idx,4);
        uncorrected_init_j = saved_final_state_init_mult_cross(uncorrected_init_idx,5);
        uncorrected_final_i = single_cross(uncorrected_final_idx,4);
        uncorrected_final_j = single_cross(uncorrected_final_idx,5);
        uncorrected_init_time = saved_final_state_init_mult_cross(uncorrected_init_idx,6);
        % uncorrected_final_time = mult_cross;
    elseif find(option_array==min(option_array)) == 4
        uncorrected_init_idx = mult_mult(2);
        uncorrected_final_idx = mult_mult(3);
        uncorrected_init_i = saved_final_state_init_mult_cross(uncorrected_init_idx,4);
        uncorrected_init_j = saved_final_state_init_mult_cross(uncorrected_init_idx,5);
        uncorrected_final_i = mult_cross(uncorrected_final_idx,4);
        uncorrected_final_j = mult_cross(uncorrected_final_idx,5);
        uncorrected_init_time = saved_final_state_init_mult_cross(uncorrected_init_idx,6);
        uncorrected_final_time = mult_cross(uncorrected_final_idx,6);
    end

    uncorrected_init_state0 = xout_lyapunov_init(uncorrected_init_j,:);
    uncorrected_final_state0 = xout_lyapunov_final(uncorrected_final_j,:);
    uncorrected_init_thrust = thrust_direction_50(:,uncorrected_init_i);
    if max(mult_cross(:,4)) > 25
        uncorrected_final_thrust = thrust_direction_50(:,uncorrected_final_i);
    else
        uncorrected_final_thrust = thrust_direction_25(:,uncorrected_final_i);
    end
    
    figure()
    scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
    hold on
    scatter(l1_pos(1), l1_pos(2), 'filled', 'green')
    plot(xout_lyapunov_init(1:uncorrected_init_j,1), xout_lyapunov_init(1:uncorrected_init_j,2), 'blue', 'LineWidth', 2)
    plot(xout_lyapunov_final(uncorrected_final_j:end,1), xout_lyapunov_final(uncorrected_final_j:end,2), 'red', 'LineWidth', 2)
    grid on
    title("Uncorrected Trajectory Option 1")
    xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
    
    if find(option_array==min(option_array)) == 1
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
        [init_output_time, uncorrected_init_xout] = ode113(fun, [0, 45], [uncorrected_init_state0, 1], options_single_cross);
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
        [final_output_time, uncorrected_final_xout] = ode113(fun, [0, -25], [uncorrected_final_state0, final_mass(i)], options_single_cross);
    elseif find(option_array==min(option_array)) == 2
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
        [init_output_time, uncorrected_init_xout] = ode113(fun, [0, 45], [uncorrected_init_state0, 1], options_single_cross);
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
        [final_output_time, uncorrected_final_xout] = ode113(fun, [0, uncorrected_final_time], [uncorrected_final_state0, final_mass(i)], options_mult_cross);
    elseif find(option_array==min(option_array)) == 3
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
        [init_output_time, uncorrected_init_xout] = ode113(fun, [0, uncorrected_init_time], [uncorrected_init_state0, 1], options_mult_cross);
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
        [final_output_time, uncorrected_final_xout] = ode113(fun, [0, -25], [uncorrected_final_state0, final_mass(i)], options_single_cross);
    elseif find(option_array==min(option_array)) == 4
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
        [init_output_time, uncorrected_init_xout] = ode113(fun, [0, uncorrected_init_time], [uncorrected_init_state0, 1], options_mult_cross);
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
        [final_output_time, uncorrected_final_xout] = ode113(fun, [0, uncorrected_final_time], [uncorrected_final_state0, final_mass(i)], options_mult_cross);
    end
    
    plot(uncorrected_init_xout(:,1), uncorrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)
    plot(uncorrected_final_xout(:,1), uncorrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
    hold off
    grid on
    legend("L2", "L1", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")

    % Correction
    % Mass burn rate
    % Non-dim thrust
    f = (T*t_star_em^2)/(l_star_em*init_mass);
    % This is the mass at the end of the each burn
    mass_1_0 = 1;
    mass_2_0 = mass_1_0;
    % mass_3_0 = saved_final_state_final_single_cross(uncorrected_final_idx, 6);
    mass_3_0 = uncorrected_final_xout(end,7);
    % mass_3_0 = mass_2_0 + mdot*abs(uncorrected_final_tout(end));
    % mass_4_0 = mass_3_0 + mdot*(tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j));
    mass_4_0 = final_mass(i,1);

    pre_moon_int = floor(length(uncorrected_init_xout)/num_pre_moon);
    post_moon_int = floor(length(uncorrected_final_xout)/num_post_moon);

    % V0 = cell(num_pre_moon+num_post_moon+2, 1);
    V_config = [0, 1];
    V1{1,1} = [xout_lyapunov_init(1,:)'; mass_1_0; tout_lyapunov_init(uncorrected_init_j)];
    V4{1,1} = [xout_lyapunov_final(uncorrected_final_j,:)'; mass_4_0; tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j)];
    V2{1,1} = [uncorrected_init_state0'; mass_2_0; uncorrected_init_thrust; init_output_time(pre_moon_int)];

    final_output_time = abs(final_output_time);

    if (final_output_time(end)-final_output_time(num_post_moon*post_moon_int)) > 0
        V3{1,1} = [uncorrected_final_xout(end,:)'; uncorrected_final_thrust; final_output_time(end)-final_output_time(num_post_moon*post_moon_int)];
        V_config = [V_config 1];
        cnt = 2;
    else
        cnt = 1;
    end
    

    for j = 1:num_pre_moon
        if j == num_pre_moon
            if (init_output_time(end)-init_output_time(j*pre_moon_int)) > 0
                V2{j+1,1} = [uncorrected_init_xout(j*pre_moon_int,:)'; uncorrected_init_thrust; init_output_time(end)-init_output_time(j*pre_moon_int)];
                V_config = [V_config 1];
            end
        else
            V2{j+1,1} = [uncorrected_init_xout(j*pre_moon_int,:)'; uncorrected_init_thrust; init_output_time((j+1)*pre_moon_int)-init_output_time(j*pre_moon_int)];
            V_config = [V_config 1];
        end
        
    end
    % cnt = 1;
    % V_config = [V_config 1];
    for k = num_post_moon:-1:1
        if k == 1
            V3{cnt,1} = [uncorrected_final_xout(k*post_moon_int,:)'; uncorrected_final_thrust; final_output_time(k*post_moon_int)];
        else
            V3{cnt,1} = [uncorrected_final_xout(k*post_moon_int,:)'; uncorrected_final_thrust; final_output_time((k)*post_moon_int)-final_output_time((k-1)*post_moon_int)];
        end
        V_config = [V_config 1];
        cnt = cnt+1;
    end
    
    
    V_config = [0, ones(1, length(V2)+length(V3)), 0];
    V0 = [V1; V2; V3; V4];
    % V4 = [uncorrected_final_xout(end,1:6)'; uncorrected_final_xout(end,7); uncorrected_final_thrust; abs(uncorrected_final_tout(end))-abs(uncorrected_final_tout(79))];
    % V5 = [uncorrected_final_xout(79,1:6)'; uncorrected_final_xout(79,7); uncorrected_final_thrust; abs(uncorrected_final_tout(79))];

    % This array determines whether a particular arc is natural or thrusting
    % 0 indicates natural arc and 1 indicates thrusting arc
    % V_config = [0, 1, 1, 0];
    % V_config = [0, 1, 1, 1, 0];

    % V0 = cell(4, 1);
    % V0{1} = V1;
    % V0{2} = V2;
    % V0{3} = V3;
    % V0{4} = V4;
    % V0{5} = V5;
    % V0{6} = V6;

    system_params = [mu, t_star_em, l_star_em, T, Isp, init_mass, f, mdot];
    % plot_modular(V0, V_config, l1_pos, l2_pos, system_params);

    % mass_1_f = mass_1_0;
    % mass_2_f = mdot*V2(end) + mass_2_0;
    % mass_3_f = 0.95;
    % % mass_3_f = mdot*V3(end) + mass_3_0;
    % mass_4_f = mass_4_0;

    x_1_des = xout_lyapunov_init(1,:)';
    % x_6_des = xout_lyapunov_final(end,:)';
    x_4_des = xout_lyapunov_final(end,:)';
    % x_2_des = uncorrected_final_xout(end,1:6)';
    % x_3_des = xout_lyapunov_final(uncorrected_final_j,:)';

    % V_soln = correction_modular(V0, V_config, system_params, x_1_des, x_6_des, l1_pos, l2_pos);

    % temp = uncorrected_final_xout(end,1:6)';
    % temp(2)= uncorrected_init_xout(end,2);

    V_soln = correction_modular(V0, V_config, system_params, x_1_des, x_4_des, l1_pos, l2_pos);
    % V_soln = correction_modular(V0, V_config, system_params, x_1_des, temp, l1_pos, l2_pos);

    %% Plot

    % plot_full_modular(V_soln, V_config, l1_pos, l2_pos, system_params);
    plot_full_comparison_modular(V0, V_soln, V_config, l1_pos, l2_pos, system_params);


end



% Compare points
% single_cross_close_pts = compare_poincare_maps_fnorm_fam(saved_final_state_init_single_cross, final_states_single);
% mult_cross_close_pts = compare_poincare_maps_fnorm_fam(saved_final_state_init_mult_cross, final_states_mult);


    

