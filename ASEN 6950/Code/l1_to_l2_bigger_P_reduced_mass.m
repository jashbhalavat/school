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

l1_lyapunov_orbits = load("V_family_L1_Lyapunov.mat").V_family;
l2_lyapunov_orbits = load("V_family_L2_Lyapunov_orbits.mat").V_family;
init_orbit_idx = 25;
final_orbit_idx = 40;
em_eq_pts = load("eq_points.mat").em_eq_pts;

l2_pos = [em_eq_pts(2,:), 0];
l1_pos = [em_eq_pts(1,:), 0];
p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure(1)
scatter(l1_pos(1), l1_pos(2), 'filled', 'black')
hold on
scatter(l2_pos(1), l2_pos(2), 'filled', 'green')
scatter(p2_pos(1), p2_pos(2), 'filled', 'cyan')
axis equal

% Set options for ode113
options_no_events = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[tout_lyapunov_init, xout_lyapunov_init] = ode113(@(t,state)CR3BP(state, mu), [0, l1_lyapunov_orbits(7,init_orbit_idx)], l1_lyapunov_orbits(1:6,init_orbit_idx), options_no_events);
plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'blue', 'LineWidth',2)

[tout_lyapunov_final, xout_lyapunov_final] = ode113(@(t,state)CR3BP(state, mu), [0, l2_lyapunov_orbits(7,final_orbit_idx)], l2_lyapunov_orbits(1:6,final_orbit_idx), options_no_events);
plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'red', 'LineWidth',2)

title("Initial and Final Orbits")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
grid on
hold off
legend("L1", "L2", "Moon", "Initial Orbit", "Final Orbit")

lyapunov_init_jacobi = jacobiConstantCR3BP(xout_lyapunov_init, mu);
lyapunov_final_jacobi = jacobiConstantCR3BP(xout_lyapunov_final, mu);

%% Transfer

num_angles = 50;
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
        
        if (xout(end,2) < 0.5) & (xout(end,2) > -0.5) 
            count = count + 1;
            saved_final_state_init_single_cross(count,:) = [xout(end,2), xout(end,4), xout(end,5), i, j, xout(end,7)];
            init_traj = plot(xout(:,1), xout(:,2), 'Color', 'blue');
        end
        % end
    end
end
title("Single Cross Two-sided Trajectories")

saved_final_state_final_single_cross = [];
count = 0;

% Single Cross
for j = 1:length(xout_lyapunov_final)
    disp("Single Cross Final Traj - " + j)
    init_state_0 = [xout_lyapunov_final(j,:), 0.985];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction(:,i), f, mdot);
        [tout, xout] = ode113(fun, [0, -10], init_state_0, options_single_cross);

        % if xout(end,1) > l2_pos(1)
        if (xout(end,2) < 0.5) & (xout(end,2) > -0.5) 
            count = count + 1;
            saved_final_state_final_single_cross(count,:) = [xout(end,2), xout(end,4), xout(end,5), i, j, xout(end,7)];
            final_traj = plot(xout(:,1), xout(:,2), 'Color', 'red');
        end
        % end
    end
end

legend([init_traj, final_traj], 'Initial Trajectory', 'Final Trajectory')
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
hold off

%% extract

B1 = saved_final_state_init_single_cross(saved_final_state_init_single_cross(:,3) > -0.1 & saved_final_state_init_single_cross(:,3) < 0.1, :);
B2 = saved_final_state_final_single_cross(saved_final_state_final_single_cross(:,3) > -0.1 & saved_final_state_final_single_cross(:,3) < 0.1, :);

% Poincare Map - other representation

figure()
scatter(B1(:,1), B1(:,2), 50, B1(:,3), 'o', 'filled')
hold on
scatter(B2(:,1), B2(:,2), 50, B2(:,3), '+')
cd = colorbar;
cd.Label.Interpreter = 'Latex';
cd.Label.String = '$$\dot{y}$$';
colormap(turbo);
grid on
xlabel("y")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
title("Single Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')
legend("Initial Trajectory", "Final Trajectory")


%% Poincare Map

figure(3)
scatter3(saved_final_state_init_single_cross(:,1), saved_final_state_init_single_cross(:,2), saved_final_state_init_single_cross(:,3), 'filled', 'blue')
hold on
scatter3(saved_final_state_final_single_cross(:,1), saved_final_state_final_single_cross(:,2), saved_final_state_final_single_cross(:,3), "filled", 'red')
xlabel("y")
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
xlabel("y")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
title("Single Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')
legend("Initial Trajectory", "Final Trajectory")

%% Poinecare map - quiver

figure()
q1 = quiver(saved_final_state_init_single_cross(:,1), saved_final_state_init_single_cross(:,2), saved_final_state_init_single_cross(:,3), saved_final_state_init_single_cross(:,6));
q1.Marker = 'o';                % Set marker style to a circle
q1.MarkerEdgeColor = 'auto';    % Match arrow border color
q1.MarkerFaceColor = 'auto';    % Fill the circle to make it a dot
q1.MarkerSize = 5;              % Adjust the size of the origin dot
hold on
q2 = quiver(saved_final_state_final_single_cross(:,1), saved_final_state_final_single_cross(:,2), saved_final_state_final_single_cross(:,3), saved_final_state_final_single_cross(:,6));
q2.Marker = 'o';                % Set marker style to a circle
q2.MarkerEdgeColor = 'auto';    % Match arrow border color
q2.MarkerFaceColor = 'auto';    % Fill the circle to make it a dot
q2.MarkerSize = 5;              % Adjust the size of the origin dot
grid on
legend("Initial Trajectory","Final Trajectory")
title("Single Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')



%% Find close points

single_cross_close_pts = compare_poincare_maps(saved_final_state_init_single_cross, saved_final_state_final_single_cross);

%% Multiple Cross

options_mult_cross = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y)multiple_cross_direction_agnostic(t,y,mu));

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
        if isempty(sol.xe) == 0
            if (sol.ye(2,end) < 0.5) && (sol.ye(2,end) > -0.5)
                seconds_init = sol.xe(end);
                saved_final_state_init_mult_cross(count,:) = [sol.ye(2,end), sol.ye(4,end), sol.ye(5,end), i, j, seconds_init, sol.ye(7,end)];
                [tout, xout] = ode113(fun, [0, seconds_init], init_state_0, options_no_events);
                init_traj = plot(xout(:,1), xout(:,2), 'Color', 'b');
            end
        end

        
    end
end

saved_final_state_final_mult_cross = [];
count = 0;

% Multiple Cross
for j = 1:length(xout_lyapunov_final)
    disp("Multiple Cross Final Traj - " + j)
    init_state_0 = [xout_lyapunov_final(j,:), 0.9];
    for i = 1:num_angles
        count = count + 1;
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, thrust_direction(:,i), f, mdot);
        sol = ode113(fun, [0, -15], init_state_0, options_mult_cross);
        if isempty(sol.xe) == 0
            if length(sol.xe) >= 2
                if (sol.ye(2,end) < 0.5) && (sol.ye(2,end) > -0.5)
                    seconds_final = sol.xe(end);
                    saved_final_state_final_mult_cross(count,:) = [sol.ye(2,end), sol.ye(4,end), sol.ye(5,end), i, j, seconds_final, sol.ye(7,end)];
                    
                    [tout, xout] = ode113(fun, [0, seconds_final], init_state_0, options_no_events);
                    final_traj = plot(xout(:,1), xout(:,2), 'Color', 'r');
                end
            end
        end
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
xlabel("y")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
zlabel("$$\dot{y}$$", 'Interpreter', 'latex')
legend("Initial Trajectory", "Final Trajectory")
title("Multi-Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')

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
xlabel("y")
ylabel("$$\dot{x}$$", 'Interpreter', 'latex')
title("Multi-Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')
legend("Initial Trajectory", "Final Trajectory")

%% Poincare map - quiver

figure()
q1 = quiver(saved_final_state_init_mult_cross(:,1), saved_final_state_init_mult_cross(:,2), saved_final_state_init_mult_cross(:,3), saved_final_state_init_mult_cross(:,7));
q1.Marker = 'o';                % Set marker style to a circle
q1.MarkerEdgeColor = 'auto';    % Match arrow border color
q1.MarkerFaceColor = 'auto';    % Fill the circle to make it a dot
q1.MarkerSize = 5;              % Adjust the size of the origin dot
hold on
q2 = quiver(saved_final_state_final_mult_cross(:,1), saved_final_state_final_mult_cross(:,2), saved_final_state_final_mult_cross(:,3), saved_final_state_final_mult_cross(:,7));
q2.Marker = 'o';                % Set marker style to a circle
q2.MarkerEdgeColor = 'auto';    % Match arrow border color
q2.MarkerFaceColor = 'auto';    % Fill the circle to make it a dot
q2.MarkerSize = 5;              % Adjust the size of the origin dot
grid on
legend("Initial Trajectory","Final Trajectory")
title("Single Cross Two-Sided Poincar\'e Map", 'Interpreter','latex')

%% Find close points

mult_cross_close_pts = compare_poincare_maps(saved_final_state_init_mult_cross, saved_final_state_final_mult_cross);

%% Plot uncorrected trajectory - option 1, single cross

% Found these nubmers by manually looking through plot:
% find(abs(saved_final_state_init_single_cross(:,1)-(-0.00693996)) < 1e-5)
% = 922
% find(abs(saved_final_state_final_single_cross(:,1)-(-0.00693925)) < 1e-5)
% = 1737
uncorrected_init_idx = single_cross_close_pts(2);
uncorrected_final_idx = single_cross_close_pts(3);
% uncorrected_init_idx = mult_cross_close_pts(2);
% uncorrected_final_idx = mult_cross_close_pts(3);
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
uncorrected_init_time = saved_final_state_init_single_cross(uncorrected_init_idx,6);
uncorrected_final_time = saved_final_state_final_single_cross(uncorrected_final_idx,6);

figure(8)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
scatter(l1_pos(1), l1_pos(2), 'filled', 'green')
plot(xout_lyapunov_init(1:uncorrected_init_j,1), xout_lyapunov_init(1:uncorrected_init_j,2), 'blue', 'LineWidth', 2)
plot(xout_lyapunov_final(uncorrected_final_j:end,1), xout_lyapunov_final(uncorrected_final_j:end,2), 'red', 'LineWidth', 2)
grid on
title("Uncorrected Trajectory Option 1")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)


fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_init_thrust, f, mdot);
[uncorrected_init_tout, uncorrected_init_xout] = ode113(fun, [0, 45], [uncorrected_init_state0, 1], options_single_cross);
x_1_f = xout(end,:)';
plot(uncorrected_init_xout(:,1), uncorrected_init_xout(:,2), 'Color', 'black', 'LineWidth', 2)

fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uncorrected_final_thrust, f, mdot);
[uncorrected_final_tout, uncorrected_final_xout] = ode113(fun, [0, -25], [uncorrected_final_state0, 0.985], options_single_cross);
plot(uncorrected_final_xout(:,1), uncorrected_final_xout(:,2), 'Color', 'magenta', 'LineWidth', 2)
hold off
grid on
legend("L2", "L1", "Initial Orbit", "Final Orbit", "Initial Transfer", "Final Transfer")

%% Correction

% Mass burn rate
% Non-dim thrust
f = (T*t_star_em^2)/(l_star_em*init_mass);
% g is assumed to be 9.80665e-3 km/s^2
mdot = -(f*l_star_em)/(Isp*9.80665e-3*t_star_em);

% This is the mass at the end of the each burn
mass_1_0 = 1;
mass_2_0 = mass_1_0;
% mass_3_0 = saved_final_state_final_single_cross(uncorrected_final_idx, 6);
mass_3_0 = uncorrected_final_xout(end,7);
% mass_3_0 = mass_2_0 + mdot*abs(uncorrected_final_tout(end));
% mass_4_0 = mass_3_0 + mdot*(tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j));
mass_4_0 = 0.985;

% V1 = [xout_lyapunov_init(1,:)'; mass_1_0; tout_lyapunov_init(uncorrected_init_j)];
% V6 = [xout_lyapunov_final(uncorrected_final_j,:)'; mass_4_0; tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j)];
% V2 = [uncorrected_init_state0'; mass_2_0; uncorrected_init_thrust; uncorrected_init_tout(74)];
% V3 = [uncorrected_init_xout(74,1:6)'; uncorrected_init_xout(74,7); uncorrected_init_thrust; uncorrected_init_tout(end)-uncorrected_init_tout(74)];
% V4 = [uncorrected_final_xout(end,1:6)'; mass_3_0; uncorrected_final_thrust; uncorrected_final_tout(76)-(uncorrected_final_tout(end))];
% V5 = [uncorrected_final_xout(76,1:6)'; uncorrected_final_xout(76,7); uncorrected_final_thrust; abs(uncorrected_final_tout(76))];

V1 = [xout_lyapunov_init(1,:)'; mass_1_0; tout_lyapunov_init(uncorrected_init_j)];
V4 = [xout_lyapunov_final(uncorrected_final_j,:)'; mass_4_0; tout_lyapunov_final(end) - tout_lyapunov_final(uncorrected_final_j)];
V2 = [uncorrected_init_state0'; mass_2_0; uncorrected_init_thrust; uncorrected_init_tout(end)];
% V3 = [uncorrected_init_xout(74,1:6)'; uncorrected_init_xout(74,7); uncorrected_init_thrust; uncorrected_init_tout(end)-uncorrected_init_tout(74)];
V3 = [uncorrected_final_xout(end,1:6)'; mass_3_0; uncorrected_final_thrust; abs(uncorrected_final_tout(end))];
% V5 = [uncorrected_final_xout(76,1:6)'; uncorrected_final_xout(76,7); uncorrected_final_thrust; abs(uncorrected_final_tout(76))];

% This array determines whether a particular arc is natural or thrusting
% 0 indicates natural arc and 1 indicates thrusting arc
% V_config = [0, 1, 1, 1, 1, 0];
V_config = [0, 1, 1, 0];

V0 = cell(4, 1);
V0{1} = V1;
V0{2} = V2;
V0{3} = V3;
V0{4} = V4;
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
