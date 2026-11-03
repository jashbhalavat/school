function output = varying_final_mass(final_nondim_mass)
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
        
    % Set options for ode113
    options_no_events = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    
    [tout_lyapunov_init, xout_lyapunov_init] = ode113(@(t,state)CR3BP(state, mu), [0, l1_lyapunov_orbits(7,init_orbit_idx)], l1_lyapunov_orbits(1:6,init_orbit_idx), options_no_events);
    
    [tout_lyapunov_final, xout_lyapunov_final] = ode113(@(t,state)CR3BP(state, mu), [0, l2_lyapunov_orbits(7,final_orbit_idx)], l2_lyapunov_orbits(1:6,final_orbit_idx), options_no_events);
    
    lyapunov_init_jacobi = jacobiConstantCR3BP(xout_lyapunov_init, mu);
    lyapunov_final_jacobi = jacobiConstantCR3BP(xout_lyapunov_final, mu);
    
    %% Transfer
    
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

    % Initial orbit, initial state
    init_state_0 = [xout_lyapunov_init(1,:), 1];

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
    
    function [value, isterminal, direction] = zero_y_hyperplane_init(t, state)
        value = state(2);
        isterminal = 0;
        direction = 0;
    end
    
    function [value, isterminal, direction] = moon_x_hyperplane_init(t, state)
        value = state(1) - (1-mu);
        isterminal = 0;
        direction = 0;
    end
    
    function [value, isterminal, direction] = zero_y_hyperplane_final(t, state)
        value = state(2);
        isterminal = 0;
        direction = 0;
    end
    
    options_single_cross = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y)single_cross_direction_agnostic(t,y,mu));

    saved_final_state_init_single_cross = load('saved_final_state_init_single_cross.mat').saved_final_state_init_single_cross;
    saved_final_state_init_mult_cross = load('saved_final_state_init_mult_cross.mat').saved_final_state_init_mult_cross;

    saved_final_state_final_single_cross = [];
    count = 0;
    
    % Single Cross
    for j = 1:length(xout_lyapunov_final)
        disp("Single Cross Final Traj - " + j)
        init_state_0 = [xout_lyapunov_final(j,:), final_nondim_mass];
        for i = 1:num_angles
            % Start from pointing straight down and rotate ccw
            fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction_25(:,i), f, mdot);
            [tout, xout] = ode113(fun, [0, -10], init_state_0, options_single_cross);
    
            % if xout(end,1) > l2_pos(1)
            % if (xout(end,2) < 0.5) & (xout(end,2) > -0.5) 
            count = count + 1;
            saved_final_state_final_single_cross(count,:) = [xout(end,2), xout(end,4), xout(end,5), i, j, xout(end,7)];
            % end
            % end
        end
    end
    
    %% Find close points
    
    single_cross_close_pts = compare_poincare_maps_fnorm(saved_final_state_init_single_cross, saved_final_state_final_single_cross);
    
    %% Multiple Cross
    
    options_mult_cross = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y)multiple_cross_direction_agnostic(t,y,mu));
    
    saved_final_state_final_mult_cross = [];
    count = 0;
    
    % Multiple Cross
    for j = 1:length(xout_lyapunov_final)
        disp("Multiple Cross Final Traj - " + j)
        init_state_0 = [xout_lyapunov_final(j,:), final_nondim_mass];
        for i = 1:num_angles
            count = count + 1;
            % Start from pointing straight down and rotate ccw
            fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction_25(:,i), f, mdot);
            sol = ode113(fun, [0, -15], init_state_0, options_mult_cross);
            if isempty(sol.xe) == 0
                if length(sol.xe) >= 2
                    % if (sol.ye(2,end) < 0.5) && (sol.ye(2,end) > -0.5)
                    seconds_final = sol.xe(end);
                    saved_final_state_final_mult_cross(count,:) = [sol.ye(2,end), sol.ye(4,end), sol.ye(5,end), i, j, seconds_final, sol.ye(7,end)];
                    % [tout, xout] = ode113(fun, [0, seconds_final], init_state_0, options_no_events);
                    % end
                end
            end
        end
    end

    % Only keep rows that are not all 0s
    saved_final_state_final_mult_cross(all(saved_final_state_final_mult_cross == 0, 2), :) = [];

    fileName = sprintf("final_nondim_mass_%d.mat", final_nondim_mass);
    save(fileName, 'saved_final_state_final_single_cross' ,'saved_final_state_final_mult_cross');
    
    %% Find close points
    
    mult_cross_close_pts = compare_poincare_maps_fnorm(saved_final_state_init_mult_cross, saved_final_state_final_mult_cross);
    
    %% Uncorrected trajectory
    
    if single_cross_close_pts(1) < mult_cross_close_pts(1)
        uncorrected_init_idx = single_cross_close_pts(2);
        uncorrected_final_idx = single_cross_close_pts(3);
        uncorrected_init_i = saved_final_state_init_single_cross(uncorrected_init_idx,4);
        uncorrected_init_j = saved_final_state_init_single_cross(uncorrected_init_idx,5);
        uncorrected_final_i = saved_final_state_final_single_cross(uncorrected_final_idx,4);
        uncorrected_final_j = saved_final_state_final_single_cross(uncorrected_final_idx,5);
        uncorrected_init_time = saved_final_state_init_single_cross(uncorrected_init_idx,6);
        uncorrected_final_time = saved_final_state_final_single_cross(uncorrected_final_idx,6);
    else
        uncorrected_init_idx = mult_cross_close_pts(2);
        uncorrected_final_idx = mult_cross_close_pts(3);
        uncorrected_init_i = saved_final_state_init_mult_cross(uncorrected_init_idx,4);
        uncorrected_init_j = saved_final_state_init_mult_cross(uncorrected_init_idx,5);
        uncorrected_final_i = saved_final_state_final_mult_cross(uncorrected_final_idx,4);
        uncorrected_final_j = saved_final_state_final_mult_cross(uncorrected_final_idx,5);
        uncorrected_init_time = saved_final_state_init_mult_cross(uncorrected_init_idx,6);
        uncorrected_final_time = saved_final_state_final_mult_cross(uncorrected_final_idx,6);
    end
    
    uncorrected_init_state0 = xout_lyapunov_init(uncorrected_init_j,:);
    uncorrected_final_state0 = xout_lyapunov_final(uncorrected_final_j,:);
    uncorrected_init_thrust = thrust_direction_50(:,uncorrected_init_i);
    uncorrected_final_thrust = thrust_direction_25(:,uncorrected_final_i);
    
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
    
    if single_cross_close_pts(1) < mult_cross_close_pts(1)
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
        [~, uncorrected_init_xout] = ode113(fun, [0, 45], [uncorrected_init_state0, 1], options_single_cross);
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
        [~, uncorrected_final_xout] = ode113(fun, [0, -25], [uncorrected_final_state0, final_nondim_mass], options_single_cross);
    else
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
        [~, uncorrected_init_xout] = ode113(fun, [0, uncorrected_init_time], [uncorrected_init_state0, 1], options_mult_cross);
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
        [~, uncorrected_final_xout] = ode113(fun, [0, uncorrected_final_time], [uncorrected_final_state0, final_nondim_mass], options_mult_cross);
    end
    
    plot(uncorrected_init_xout(:,1), uncorrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)
    plot(uncorrected_final_xout(:,1), uncorrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
    hold off
    grid on
    legend("L2", "L1", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")

    output = [single_cross_close_pts(1), mult_cross_close_pts(1), uncorrected_init_xout(end,7), uncorrected_final_xout(end,7)];
    end