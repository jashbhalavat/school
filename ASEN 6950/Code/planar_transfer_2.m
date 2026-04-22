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

% Non-dim thrust
f = (T*t_star_em^2)/(l_star_em*init_mass);

% g is assumed to be 9.80665e-3 km/s^2
mdot = -(f*l_star_em)/(Isp*9.80665e-3*t_star_em);

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
init_state_0 = [xout_lyapunov_init(1,:), 1];

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
    init_state_0 = [xout_lyapunov_init(j,:), 1];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction(:,i), f, mdot);
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
    init_state_0 = [xout_lyapunov_final(j,:), 1];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction(:,i), f, mdot);
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
    init_state_0 = [xout_lyapunov_init(j,:), 1];
    for i = 1:num_angles
        count = count + 1;
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction(:,i), f, mdot);
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
    init_state_0 = [xout_lyapunov_final(j,:), 1];
    for i = 1:num_angles
        count = count + 1;

        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction(:,i), f, mdot);
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

mult_cross_close_pts = compare_poincare_maps_beyond_moon(saved_final_state_init_mult_cross, saved_final_state_final_mult_cross, p2_pos(1));
% mult_cross_close_pts = compare_poincare_maps(saved_final_state_init_mult_cross, saved_final_state_final_mult_cross);
        

%% Plot uncorrected trajectory - option 2, single cross

% These numbers were found by the following commands:
% [~, idx] = min(abs(saved_final_state_final_mult_cross(:,1) - -3.474))
% [~, idx] = min(abs(saved_final_state_init_mult_cross(:,1) - -3.474))

uncorrected_init_idx = 246;
uncorrected_final_idx = 3418;
uncorrected_init_i = saved_final_state_init_mult_cross(uncorrected_init_idx,4);
uncorrected_init_j = saved_final_state_init_mult_cross(uncorrected_init_idx,5);
uncorrected_init_time = saved_final_state_init_mult_cross(uncorrected_init_idx,6);
uncorrected_final_i = saved_final_state_final_mult_cross(uncorrected_final_idx,4);
uncorrected_final_j = saved_final_state_final_mult_cross(uncorrected_final_idx,5);
uncorrected_final_time = saved_final_state_final_mult_cross(uncorrected_final_idx,6);
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
title("Uncorrected Trajectory Option 2")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)


fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
[uncorrected_init_tout, uncorrected_init_xout] = ode113(fun, [0, uncorrected_init_time], [uncorrected_init_state0, 1], options_no_events);
x_1_f = xout(end,:)';
plot(uncorrected_init_xout(:,1), uncorrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
[uncorrected_final_tout, uncorrected_final_xout] = ode113(fun, [0, uncorrected_final_time], [uncorrected_final_state0, 1], options_no_events);
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


% This is the mass at the end of the each burn
mass_1_0 = 1;
mass_2_0 = mass_1_0;
mass_3_0 = mass_2_0 + mdot*abs(uncorrected_final_tout(end));
mass_4_0 = mass_3_0 + mdot*(tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j));

V1 = [xout_lyapunov_init(1,:)'; mass_1_0; tout_lyapunov_init(uncorrected_init_j)];
V4 = [xout_lyapunov_final(uncorrected_final_j,:)'; mass_4_0; tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j-1)];
V2 = [uncorrected_init_state0'; mass_2_0; uncorrected_init_thrust; uncorrected_init_tout(end)];
V3 = [uncorrected_final_xout(end,1:6)'; mass_3_0; uncorrected_final_thrust; abs(uncorrected_final_tout(end-1))];

V0 = [V1; V2; V3; V4];

mass_1_f = mass_1_0;
mass_2_f = mdot*V0(19) + mass_2_0;
mass_3_f = mdot*V0(30) + mass_3_0;
mass_4_f = mass_4_0;

x_1_des = xout_lyapunov_init(1,:)';
x_4_des = xout_lyapunov_final(end,:)';
x_2_des = uncorrected_final_xout(end,1:6)';
x_3_des = xout_lyapunov_final(uncorrected_final_j,:)';

system_params = [mu, t_star_em, l_star_em, T, Isp, init_mass, f, mdot];

V_soln = correction(V0, system_params, x_1_des, x_2_des, x_3_des, x_4_des);

%% Plot corrected trajectory - first two arcs

figure(9)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
% plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'blue', 'LineWidth', 2)
% plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'red', 'LineWidth', 2)
grid on
title("Uncorrected Trajectory Option 1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)

[corrected_tout_1, corrected_xout_1] = ode113(@(t,state)CR3BP(state, mu), [0, V_soln(7)], V_soln(1:6), options_no_events);
plot(corrected_xout_1(:,1), corrected_xout_1(:,2), 'Color', 'Blue', 'LineWidth', 2)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, V_soln(16:18), f, mdot);
[corrected_tout_2, corrected_xout_2] = ode113(fun, [0, V_soln(15)], [V_soln(8:13); V_soln(14)], options_no_events);
plot(corrected_xout_2(:,1), corrected_xout_2(:,2), 'Color', 'Black', 'LineWidth', 2)

% corrected_final_mass = init_mass + mdot*V_soln(7) - mdot*V_soln(14);
% 
% fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, init_mass);
% [corrected_final_tout, corrected_final_xout] = ode113(fun, [0, V_soln(14)], [V_soln(8:13); corrected_final_mass], options_no_events);
% plot(corrected_final_xout(:,1), corrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
hold off
grid on
legend("L2", "Initial Orbit", "Initial Transfer")

figure(10)
plot(corrected_tout_1, ones(length(corrected_tout_1),1)*init_mass, 'Color', 'Blue', 'LineWidth', 2)
hold on
grid on
xlabel("Non-dimensional Time")
ylabel("Spacecraft Mass [kg]")
title("Spacecraft Mass over Time")
plot(corrected_tout_1(end) + corrected_tout_2, corrected_xout_2(:,7), 'Color', 'Black', 'LineWidth', 2)
legend("Initial Orbit", "Initial Transfer")
hold off

figure(11)
plot(corrected_tout_1, zeros(length(corrected_tout_1),1),'x', 'Color', 'Blue', 'LineWidth', 2)
plot(corrected_tout_1, zeros(length(corrected_tout_1),1),'+', 'Color', 'Blue', 'LineWidth', 2)
plot(corrected_tout_1, zeros(length(corrected_tout_1),1),'*', 'Color', 'Blue', 'LineWidth', 2)
hold on
grid on
xlabel("Non-dimensional Time")
ylabel("Thrust Vector")
title("Spacecraft Thrust Vector over Time")
plot(corrected_tout_1(end) + corrected_tout_2, ones(length(corrected_tout_2),3)*V_soln(9),'x', 'Color', 'Black', 'LineWidth', 2)
plot(corrected_tout_1(end) + corrected_tout_2, ones(length(corrected_tout_2),3)*V_soln(10),'--', 'Color', 'Black', 'LineWidth', 2)
plot(corrected_tout_1(end) + corrected_tout_2, ones(length(corrected_tout_2),3)*V_soln(11), '-.', 'Color', 'Black', 'LineWidth', 2)
legend('$$\hat{u}_{1,x}$$', '$$\hat{u}_{1,y}$$', '$$\hat{u}_{1,z}$$', '$$\hat{u}_{2,x}$$', '$$\hat{u}_{2,y}$$', '$$\hat{u}_{2,z}$$','Interpreter','Latex')
hold off

%% Plot corrected trajectory - first four arcs

figure(12)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
plot(xout_lyapunov_init(1:uncorrected_init_j,1), xout_lyapunov_init(1:uncorrected_init_j,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 4)
plot(xout_lyapunov_final(uncorrected_final_j:end,1), xout_lyapunov_final(uncorrected_final_j:end,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
grid on
title("Uncorrected Trajectory Option 1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
plot(uncorrected_init_xout(:,1), uncorrected_init_xout(:,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
plot(uncorrected_final_xout(:,1), uncorrected_final_xout(:,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
grid on
legend("L2", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")

title("Corrected Trajectory Option 1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)

[corrected_tout_1, corrected_xout_1] = ode113(@(t,state)CR3BP(state, mu), [0, V_soln(8)], V_soln(1:6), options_no_events);
plot(corrected_xout_1(:,1), corrected_xout_1(:,2), '--', 'Color', 'Blue', 'LineWidth', 2)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, V_soln(16:18), f, mdot);
[corrected_tout_2, corrected_xout_2] = ode113(fun, [0, V_soln(19)], V_soln(9:15), options_no_events);
plot(corrected_xout_2(:,1), corrected_xout_2(:,2), '--', 'Color', 'Black', 'LineWidth', 2)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, V_soln(27:29), f, mdot);
[corrected_tout_3, corrected_xout_3] = ode113(fun, [0, V_soln(30)], V_soln(20:26), options_no_events);
plot(corrected_xout_3(:,1), corrected_xout_3(:,2), '--', 'Color', 'Magenta', 'LineWidth', 2)

[corrected_tout_4, corrected_xout_4] = ode113(@(t,state)CR3BP(state, mu), [0, V_soln(38)], V_soln(31:36), options_no_events);
plot(corrected_xout_4(:,1), corrected_xout_4(:,2), '--', 'Color', 'Red', 'LineWidth', 2)

hold off
legend("L2", "Uncorrected Trajectory", "Initial Orbit", "Initial Transfer", "Final Transfer", "Final Orbit")

%% Plot corrected mass

figure(13)
plot([0, V0(8)], [mass_1_0, mass_1_f], 'Color', [0.5 0.5 0.5], 'LineWidth', 4)
hold on
plot([V0(8), V0(8)+V0(19)], [mass_2_0, mass_2_f], 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
plot([V0(8)+V0(19), V0(8)+V0(19)+V0(30)], [mass_3_0, mass_3_f], 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
plot([V0(8)+V0(19)+V0(30), V0(8)+V0(19)+V0(30)+V0(38)], [mass_4_0, mass_4_f], 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
plot(corrected_tout_1, ones(length(corrected_tout_1),1), '--', 'Color', 'Blue', 'LineWidth', 2)
grid on
xlabel("Non-dimensional Time")
ylabel("Spacecraft Non-Dimensional Mass")
title("Spacecraft Mass over Time")
plot(corrected_tout_1(end) + corrected_tout_2, corrected_xout_2(:,7), '--', 'Color', 'Black', 'LineWidth', 2)
plot(corrected_tout_1(end) + corrected_tout_2(end) + corrected_tout_3, corrected_xout_3(:,7), '--', 'Color', 'Magenta', 'LineWidth', 2)
plot(corrected_tout_1(end) + corrected_tout_2(end) + corrected_tout_3(end) + corrected_tout_4, ones(length(corrected_tout_4),1)*V_soln(37), '--', 'Color', 'Red', 'LineWidth', 2)
legend("Uncorrected Mass", "Initial Orbit", "Initial Transfer", "Final Transfer", "Final Orbit")
hold off

%% Plot Jacobi Constant

figure(14)
c_uncorrected_1 = jacobiConstantCR3BP(xout_lyapunov_init, mu);
c_uncorrected_2 = jacobiConstantCR3BP(uncorrected_init_xout, mu);
c_uncorrected_3 = jacobiConstantCR3BP(uncorrected_final_xout, mu);
c_uncorrected_4 = jacobiConstantCR3BP(xout_lyapunov_final, mu);

uncorrected_time = [tout_lyapunov_init; tout_lyapunov_init(end)+uncorrected_init_tout;
                    tout_lyapunov_init(end)+uncorrected_init_tout(end)+uncorrected_final_tout;
                    tout_lyapunov_init(end)+uncorrected_init_tout(end)+uncorrected_final_tout(end)+tout_lyapunov_final];

c_corrected_1 = jacobiConstantCR3BP(corrected_xout_1, mu);
c_corrected_2 = jacobiConstantCR3BP(corrected_xout_2, mu);
c_corrected_3 = jacobiConstantCR3BP(corrected_xout_3, mu);
c_corrected_4 = jacobiConstantCR3BP(corrected_xout_4, mu);

plot(tout_lyapunov_init(1:uncorrected_init_j), c_uncorrected_1(1:uncorrected_init_j), 'Color', [0.5 0.5 0.5], 'LineWidth', 4)
hold on
plot(tout_lyapunov_init(uncorrected_init_j)+uncorrected_init_tout, c_uncorrected_2, 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
plot(tout_lyapunov_init(uncorrected_init_j)+uncorrected_init_tout(end)+uncorrected_final_tout, c_uncorrected_3, 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
plot(tout_lyapunov_init(uncorrected_init_j)+uncorrected_init_tout(end)+uncorrected_final_tout(end)+tout_lyapunov_final(uncorrected_final_j:end), c_uncorrected_4(uncorrected_final_j:end), 'Color', [0.5 0.5 0.5], 'LineWidth', 4, 'HandleVisibility','off')
plot(corrected_tout_1, c_corrected_1, '--', 'Color', 'Blue', 'LineWidth', 2)
plot(corrected_tout_1(end) + corrected_tout_2, c_corrected_2, '--', 'Color', 'Black', 'LineWidth', 2)
plot(corrected_tout_1(end) + corrected_tout_2(end) + corrected_tout_3, c_corrected_3, '--', 'Color', 'Magenta', 'LineWidth', 2)
plot(corrected_tout_1(end) + corrected_tout_2(end) + corrected_tout_3(end) + corrected_tout_4, c_corrected_4, '--', 'Color', 'Red', 'LineWidth', 2)
grid on
xlabel("Non-dimensional Time")
ylabel("Spacecraft Non-Dimensional Mass")
title("Spacecraft Jacobi Constant over Time")
legend("Uncorrected Mass", "Initial Orbit", "Initial Transfer", "Final Transfer", "Final Orbit")
hold off



%% Mdot measured

delta_time = corrected_tout_1(end) + corrected_tout_2(end) + corrected_tout_3(end) - corrected_tout_1(end);
mass_diff = V_soln(36) - 1;
mdot_measured = mass_diff/delta_time;


%% Plot corrected trajectory - option 1, single cross

figure(9)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
% plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'blue', 'LineWidth', 2)
% plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'red', 'LineWidth', 2)
grid on
title("Uncorrected Trajectory Option 1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
[corrected_init_tout, corrected_init_xout] = ode113(fun, [0, V_soln(7)], [V_soln(1:6); 1], options_no_events);
plot(corrected_init_xout(:,1), corrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)

corrected_final_mass = init_mass + mdot*V_soln(7) - mdot*V_soln(14);

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
[corrected_final_tout, corrected_final_xout] = ode113(fun, [0, V_soln(14)], [V_soln(8:13); corrected_final_mass], options_no_events);
plot(corrected_final_xout(:,1), corrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
hold off
grid on
legend("L2", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")



%% Jacobi constant


