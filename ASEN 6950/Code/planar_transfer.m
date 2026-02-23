clear; clc; close all;

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

lyapunov_orbits = load("V_family_L2_Lyapunov_orbits.mat").V_family;
init_orbit_idx = 5;
final_orbit_idx = 10;
em_eq_pts = load("eq_points.mat").em_eq_pts;

l2_pos = [em_eq_pts(2,:), 0];
l1_pos = [em_eq_pts(1,:), 0];
p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure(1)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on

% Set options for ode113
options_no_events = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[tout_lyapunov_init, xout_lyapunov_init] = ode113(@(t,state)CR3BP(state, mu), [0, lyapunov_orbits(7,init_orbit_idx)], lyapunov_orbits(1:6,init_orbit_idx), options_no_events);
plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'blue', 'LineWidth',2)

[tout_lyapunov_final, xout_lyapunov_final] = ode113(@(t,state)CR3BP(state, mu), [0, lyapunov_orbits(7,final_orbit_idx)], lyapunov_orbits(1:6,final_orbit_idx), options_no_events);
plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'red', 'LineWidth',2)

title("Initial and Final Orbits")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
grid on
hold off
legend("L2", "Initial Orbit", "Final Orbit")

%% Transfer

num_angles = 25;
angles = linspace(0, 2*pi, num_angles);

% Straight down vector
neg_y_vec = [0, -1, 0]';

for i = 1:num_angles
    dcm = R3(angles(i));
    thrust_direction(:,i) = dcm * neg_y_vec;
end

% Initial orbit, initial state
init_state_0 = [xout_lyapunov_init(1,:), init_mass];

function status = countNegCrossings(t,state,flag)
    persistent count lastSign
    status = 0;

    if strcmp(flag,'init')
        count = 0;
        lastSign = sign(state(2,1));
    elseif isempty(flag)
        s = sign(state(2,end));
        if (lastSign > 0 && s < 0) || (lastSign < 0 && s > 0)
            count = count + 1;
            if count == 2
                status = 1;  % STOP
            end
        end
        lastSign = s;
    end
end

function status = countPosCrossings(t,state,flag)
    persistent count lastSign
    status = 0;

    if strcmp(flag,'init')
        count = 0;
        lastSign = sign(state(2,1));
    elseif isempty(flag)
        s = sign(state(2,end));
        if (lastSign < 0 && s > 0) || (lastSign > 0 && s < 0)
            count = count + 1;
            if count == 2
                status = 1;  % STOP
            end
        end
        lastSign = s;
    end
end

function [value, isterminal, direction] = single_cross_direction_agnostic(t, y)
    value = y(2);    % Detect when y(2) crosses 0
    direction = 0;   % Trigger for both increasing and decreasing
    isterminal = 1;  % Stop after one event
end

function [value, isterminal, direction] = multiple_cross_direction_agnostic(t, y)
    value = y(2);    % Detect when y(2) crosses 0
    direction = 0;   % Trigger for both increasing and decreasing
    isterminal = 0;  % Don't stop at event trigger
end

function [value, isterminal, direction] = zero_y_hyperplane_init(t, state)
    value = state(2);
    isterminal = 0;
    direction = 0;
    % persistent eventCount
    % if isempty(crossedOnce)
    %     crossedOnce = false;
    % end
    % 
    % if ~crossedOnce
    %     % First zero crossing (non-terminal)
    %     value      = [state(2); 1];
    %     isterminal = [0; 0];
    % else
    %     % Second zero crossing (terminal)
    %     value      = [1; state(2)];
    %     isterminal = [0; 1];
    % end
    % 
    % direction = [-1; -1];
    % 
    % % Detect completed crossing
    % if state(2) == 0
    %     crossedOnce = true;
    % end

    % if isempty(eventCount)
    %     eventCount = 0;
    % end
    % 
    % value = state(2);        % detect y(2) = 0
    % if value == 0
    %     eventCount = eventCount + 1;
    % end
    % direction = 0;      % ONLY positive → negative crossings
    % 
    % eventCount = eventCount + 1;
    % 
    % if eventCount >= 2
    %     isterminal = 1;  % stop after 2nd crossing
    % else
    %     isterminal = 0;  % keep going
    % end

    % x_cross_count = x_cross_count + 1;
    % 
    % if x_cross_count < 2
    %     isterminal = 0;
    % else
    %     isterminal = 1;
    % end

    % persistent init_count;
    % if isempty(init_count)
    %     init_count = 0;
    % end


    % value = state(2);          % Stop when y = 0
    % isterminal = 1;     % 1 = stop integration
    % direction = -1;     % Direction agnostic

    % if value == 0
    %     init_count = init_count + 1;
    % end
    % if init_count >= 2
    %     isterminal = 1;
    % else
    %     isterminal = 0;
    % end
end

function [value, isterminal, direction] = zero_y_hyperplane_final(t, state)
    % persistent count
    % if isempty(count)
    %     count = 0;
    % end
    % 
    % value = state(2);        % detect y(2) = 0
    % direction = 1;      % ONLY positive → negative crossings
    % 
    % count = count + 1;
    % 
    % if count >= 2
    %     isterminal = 1;  % stop after 2nd crossing
    % else
    %     isterminal = 0;  % keep going
    % end
    value = state(2);
    isterminal = 0;
    direction = 0;
end

options_single_cross = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @single_cross_direction_agnostic);

saved_final_state_init_single_cross = [];
count = 0;

% Single Cross
figure(2)
plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'blue', 'LineWidth',2)
hold on
grid on
plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'red', 'LineWidth',2)


for j = 1:length(xout_lyapunov_init)
    disp("Single Cross Init Traj - " + j)
    init_state_0 = [xout_lyapunov_init(j,:), init_mass];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, thrust_direction(:,i), t_star_em, l_star_em, init_mass);
        [tout, xout] = ode113(fun, [0, 10], init_state_0, options_single_cross);
        
        % if xout(end,1) > l2_pos(1)
        count = count + 1;
        saved_final_state_init_single_cross(count,:) = [xout(end,1), xout(end,4), xout(end,5), i, j];
        init_traj = plot(xout(:,1), xout(:,2), 'Color', 'blue');
        % end
    end
end
title("Single Cross Two-sided Trajectories")

saved_final_state_final_single_cross = [];
count = 0;

% Single Cross
for j = 1:length(xout_lyapunov_final)
    disp("Single Cross Final Traj - " + j)
    init_state_0 = [xout_lyapunov_final(j,:), init_mass];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, thrust_direction(:,i), t_star_em, l_star_em, init_mass);
        [tout, xout] = ode113(fun, [0, -25], init_state_0, options_single_cross);

        % if xout(end,1) > l2_pos(1)
        count = count + 1;
        saved_final_state_final_single_cross(count,:) = [xout(end,1), xout(end,4), xout(end,5), i, j];
        final_traj = plot(xout(:,1), xout(:,2), 'Color', 'red');
        % end
    end
end

legend([init_traj, final_traj], 'Initial Trajectory', 'Final Trajectory')
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
hold off

%% Poincare Map

figure(3)
scatter3(saved_final_state_init_single_cross(:,1), saved_final_state_init_single_cross(:,2), saved_final_state_init_single_cross(:,3), 'filled', 'blue')
hold on
scatter3(saved_final_state_final_single_cross(:,1), saved_final_state_final_single_cross(:,2), saved_final_state_final_single_cross(:,3), "filled", 'red')
xlabel("x")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
zlabel("$$\dot{y}$$", 'Interpreter', 'latex')
legend("Initial Trajectory", "Final Trajectory")
title("Single Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')

%% Poincare Map - other representation

figure(4)
scatter(saved_final_state_init_single_cross(:,1), saved_final_state_init_single_cross(:,2), 50, saved_final_state_init_single_cross(:,3), 'o', 'filled')
hold on
scatter(saved_final_state_final_single_cross(:,1), saved_final_state_final_single_cross(:,2), 50, saved_final_state_final_single_cross(:,3), '+')
cd = colorbar;
cd.Label.Interpreter = 'Latex';
cd.Label.String = '$$\dot{y}$$';
colormap(turbo);
grid on
xlabel("x")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
title("Single Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')
legend("Initial Trajectory", "Final Trajectory")

%% Find close points

single_cross_close_pts = compare_poincare_maps(saved_final_state_init_single_cross, saved_final_state_final_single_cross);


%% Multiple Cross

options_mult_cross = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @multiple_cross_direction_agnostic);

saved_final_state_init_mult_cross = [];
count = 0;

figure(5)
plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'blue', 'LineWidth',2)
hold on
grid on
plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'red', 'LineWidth',2)

% Multiple Cross
for j = 1:length(xout_lyapunov_init)
    disp("Multiple cross Init Traj - " + j)
    init_state_0 = [xout_lyapunov_init(j,:), init_mass];
    for i = 1:num_angles
        count = count + 1;
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, thrust_direction(:,i), t_star_em, l_star_em, init_mass);
        sol = ode113(fun, [0, 15], init_state_0, options_mult_cross);
        if init_state_0(2) == 0
            seconds_init = sol.xe(3);
            saved_final_state_init_mult_cross(count,:) = [sol.ye(1,3), xout(4,3), xout(5,3), i, j, seconds_init];
        else
            seconds_init = sol.xe(2);
            saved_final_state_init_mult_cross(count,:) = [sol.ye(1,2), xout(4,2), xout(5,2), i, j, seconds_init];
        end

        [tout, xout] = ode113(fun, [0, seconds_init], init_state_0, options_no_events);
        init_traj = plot(xout(:,1), xout(:,2), 'Color', 'b');
    end
end

saved_final_state_final_mult_cross = [];
count = 0;

% Multiple Cross
for j = 1:length(xout_lyapunov_final)
    disp("Multiple Cross Final Traj - " + j)
    init_state_0 = [xout_lyapunov_final(j,:), init_mass];
    for i = 1:num_angles
        count = count + 1;

        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, thrust_direction(:,i), t_star_em, l_star_em, init_mass);
        sol = ode113(fun, [0, -15], init_state_0, options_mult_cross);
        if init_state_0(2) == 0
            seconds_final = sol.xe(3);
            saved_final_state_final_mult_cross(count,:) = [sol.ye(1,3), xout(4,3), xout(5,3), i, j, seconds_final];
        else
            seconds_final = sol.xe(2);
            saved_final_state_final_mult_cross(count,:) = [sol.ye(1,2), xout(4,2), xout(5,2), i, j, seconds_final];
        end

        [tout, xout] = ode113(fun, [0, seconds_final], init_state_0, options_no_events);
        final_traj = plot(xout(:,1), xout(:,2), 'Color', 'r');
    end
end

hold off
legend([init_traj, final_traj], "Initial Trajectory", "Final Trajectory")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
title("2-Cross Two-sided Trajectories")


%% Mult-Cross Poincare Map

figure(6)
scatter3(saved_final_state_init_mult_cross(:,1), saved_final_state_init_mult_cross(:,2), saved_final_state_init_mult_cross(:,3), 'filled', 'blue')
hold on
scatter3(saved_final_state_final_mult_cross(:,1), saved_final_state_final_mult_cross(:,2), saved_final_state_final_mult_cross(:,3), "filled", 'red')
xlabel("x")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
zlabel("$$\dot{y}$$", 'Interpreter', 'latex')
legend("Initial Trajectory", "Final Trajectory")
title("2-Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')

%% Mult-Cross Poincare Map different representation

figure(7)
scatter(saved_final_state_init_mult_cross(:,1), saved_final_state_init_mult_cross(:,2), 50, saved_final_state_init_mult_cross(:,3), 'o', 'filled')
hold on
scatter(saved_final_state_final_mult_cross(:,1), saved_final_state_final_mult_cross(:,2), 50, saved_final_state_final_mult_cross(:,3), '+')
cd = colorbar;
cd.Label.Interpreter = 'Latex';
cd.Label.String = '$$\dot{y}$$';
colormap(turbo);
grid on
xlabel("x")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
title("2-Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')
legend("Initial Trajectory", "Final Trajectory")

%% Find close points

mult_cross_close_pts = compare_poincare_maps(saved_final_state_init_mult_cross, saved_final_state_final_mult_cross);
        

%% Plot uncorrected trajectory - option 1, single cross

uncorrected_init_idx = single_cross_close_pts(2);
uncorrected_final_idx = single_cross_close_pts(3);
uncorrected_init_i = saved_final_state_init_single_cross(uncorrected_init_idx,4);
uncorrected_init_j = saved_final_state_init_single_cross(uncorrected_init_idx,5);
% uncorrected_init_seconds = saved_final_state_init_mult_cross(uncorrected_init_idx,6);
uncorrected_final_i = saved_final_state_final_single_cross(uncorrected_final_idx,4);
uncorrected_final_j = saved_final_state_final_single_cross(uncorrected_final_idx,5);
% uncorrected_final_seconds = saved_final_state_final_mult_cross(uncorrected_final_idx,6);
uncorrected_init_state0 = xout_lyapunov_init(uncorrected_init_j,:);
uncorrected_final_state0 = xout_lyapunov_final(uncorrected_final_j,:);
uncorrected_init_thrust = thrust_direction(:,uncorrected_init_i);
uncorrected_final_thrust = thrust_direction(:,uncorrected_final_i);

figure(8)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
plot(xout_lyapunov_init(1:uncorrected_init_j,1), xout_lyapunov_init(1:uncorrected_init_j,2), 'blue', 'LineWidth', 2)
plot(xout_lyapunov_final(uncorrected_final_j:end,1), xout_lyapunov_final(uncorrected_final_j:end,2), 'red', 'LineWidth', 2)
grid on
title("Uncorrected Trajectory Option 1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)


fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_init_thrust, t_star_em, l_star_em, init_mass);
[uncorrected_init_tout, uncorrected_init_xout] = ode113(fun, [0, 25], [uncorrected_init_state0, init_mass], options_single_cross);
x_1_f = xout(end,:)';
plot(uncorrected_init_xout(:,1), uncorrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, init_mass);
[uncorrected_final_tout, uncorrected_final_xout] = ode113(fun, [0, -25], [uncorrected_final_state0, init_mass], options_single_cross);
plot(uncorrected_final_xout(:,1), uncorrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
hold off
grid on
legend("L2", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")

%% Correction

% Mass burn rate
% Non-dim thrust
f = (T*t_star_em^2)/(l_star_em*init_mass);
% g is assumed to be 9.80665e-3 km/s^2
mdot = -(f*l_star_em)/(Isp*9.80665e-3*t_star_em);

% This is the mass at the end of the burn, so it should have mass removed
% from both the init and final trajectories
uncorrected_final_mass = init_mass + mdot*uncorrected_init_tout(end) - mdot*uncorrected_final_tout(end);

V0 = [xout_lyapunov_init(1,:)'; tout_lyapunov_init(uncorrected_init_j);
    uncorrected_init_state0'; uncorrected_init_tout(end);
    uncorrected_final_state0'; uncorrected_final_tout(end);
    xout_lyapunov_final(uncorrected_final_j,:)'; tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j)];

desired_stage_1 = xout_lyapunov_init(1,:)';
desired_stage_4 = xout_lyapunov_final(end,:)';

% 
% x_1_f = uncorrected_init_xout(end,:)';
% 
% e1 = x_1_f - [uncorrected_final_state0'; uncorrected_final_mass];
% F = e1;

% V_soln = correction(V0, mu, true);
V_soln = correction(V0, mu, T, Isp, uncorrected_init_thrust, uncorrected_final_thrust, t_star_em, l_star_em, init_mass, uncorrected_final_mass, desired_stage_1, desired_stage_4);

%% Plot corrected trajectory - option 1, single cross

figure(9)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'blue', 'LineWidth', 2)
plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'red', 'LineWidth', 2)
grid on
title("Uncorrected Trajectory Option 1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_init_thrust, t_star_em, l_star_em, init_mass);
[corrected_init_tout, corrected_init_xout] = ode113(fun, [0, V_soln(7)], [V_soln(1:6); init_mass], options_no_events);
plot(corrected_init_xout(:,1), corrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)

corrected_final_mass = init_mass + mdot*V_soln(7) - mdot*V_soln(14);

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, init_mass);
[corrected_final_tout, corrected_final_xout] = ode113(fun, [0, V_soln(14)], [V_soln(8:13); corrected_final_mass], options_no_events);
plot(corrected_final_xout(:,1), corrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
hold off
grid on
legend("L2", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")



%% Jacobi constant

function C = jacobi_constant(x, mu)
    C = u_star_times_2(x(1), x(2), x(3), mu) - x(4)^2 - x(5)^2 - x(6)^2;
end

