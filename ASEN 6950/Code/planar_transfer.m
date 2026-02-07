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
m_star_em = mass_earth + mass_moon;
l_star_em = a_moon;
t_star_em = sqrt(l_star_em^3/(G * m_star_em));

mu = mass_ratio_em;

T = 15e-3; % Nm
Isp = 1320; % sec
% T_star = T/(mu*l_star_em); % non-dim
% T_star = T*t_star_em^2/l_star_em;
T_star = T;

lyapunov_orbits = load("V_family_L2_Lyapunov_orbits.mat").V_family;
em_eq_pts = load("eq_points.mat").em_eq_pts;

l2_pos = [em_eq_pts(2,:), 0];
l1_pos = [em_eq_pts(1,:), 0];
p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure(1)
scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
hold on

% Set options for ode113
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[tout_lyapunov_init, xout_lyapunov_init] = ode113(@(t,state)CR3BP(state, mu), [0, lyapunov_orbits(7,5)], lyapunov_orbits(1:6,5), options);
plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'LineWidth',2)

[tout_lyapunov_final, xout_lyapunov_final] = ode113(@(t,state)CR3BP(state, mu), [0, lyapunov_orbits(7,10)], lyapunov_orbits(1:6,10), options);
plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'LineWidth',2)

% fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, lyapunov_orbits(1:6,10));
% [tout_transfer, xout_transfer] = ode113(fun, [0, 6], [lyapunov_orbits(1:6,5); 10], options);
% plot(xout_transfer(:,1), xout_transfer(:,2), xout_transfer(:,3), 'LineWidth',2)

% legend("L2", "Initial Lyapunov Orbit", "Final Lyapunov Orbit")
title("Initial and Final Orbits")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
% zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)

%% Transfer

% Initial orbit, initial state
init_mass = 10;
init_state_0 = [xout_lyapunov_init(1,:), init_mass];

lyapunov_final_x_max = max(xout_lyapunov_final(:,1));
lyapunov_final_x_min = min(xout_lyapunov_final(:,1));

% global x_cross_count
% x_cross_count = 0;

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

% function [value, isterminal, direction] = myEvent(t, y)
%     value = y(2);    % Detect when y(2) crosses 0
%     direction = 0;         % Trigger for both increasing and decreasing
%     isterminal = 0;
% end

function [value, isterminal, direction] = single_cross_direction_agnostic(t, y)
    value = y(2);    % Detect when y(2) crosses 0
    direction = 0;         % Trigger for both increasing and decreasing
    isterminal = 1;
end

function [value, isterminal, direction] = multiple_cross_direction_agnostic(t, y)
    value = y(2);    % Detect when y(2) crosses 0
    direction = 0;         % Trigger for both increasing and decreasing
    isterminal = 0;
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

num_angles = 25;
angles = linspace(0, 2*pi, num_angles);

% Straight down vector
neg_y_vec = [0, -1, 0]';

for i = 1:num_angles
    dcm = R3(angles(i));
    thrust_direction(:,i) = dcm * neg_y_vec;
end

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @single_cross_direction_agnostic);

saved_final_state_init_single_cross = [];
count = 0;

% Single Cross
for j = 1:length(xout_lyapunov_init)
% for j = 1:1
    disp(j)
    init_state_0 = [xout_lyapunov_init(j,:), init_mass];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction(:,i));
        [tout, xout] = ode113(fun, [0, 25], init_state_0, options);

        count = count + 1;
        saved_final_state_init_single_cross(count,:) = [xout(end,1), xout(end,4), xout(end,5)];
        plot(xout(:,1), xout(:,2), 'Color', 'black')
    end
end

saved_final_state_final_single_cross = [];
count = 0;

% Single Cross
for j = 1:length(xout_lyapunov_final)
% for j = 1:1
    disp(j)
    init_state_0 = [xout_lyapunov_final(j,:), init_mass];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction(:,i));
        [tout, xout] = ode113(fun, [0, -25], init_state_0, options);

        count = count + 1;
        saved_final_state_final_single_cross(count,:) = [xout(end,1), xout(end,4), xout(end,5)];
        plot(xout(:,1), xout(:,2), 'Color', 'magenta')
    end
end

legend("L2", "Initial Lyapunov Orbit", "Final Lyapunov Orbit")

%%

figure(2)
scatter3(saved_final_state_init_single_cross(:,1), saved_final_state_init_single_cross(:,2), saved_final_state_init_single_cross(:,3), 'filled', 'black')
hold on
scatter3(saved_final_state_final_single_cross(:,1), saved_final_state_final_single_cross(:,2), saved_final_state_final_single_cross(:,3), "filled", 'magenta')
xlabel("x")
ylabel("xdot")
zlabel("ydot")
legend("Initial Orbit", "Final Orbit")
title("Single Cross Poincar\'e Map", 'Interpreter','latex')


%% Multiple Cross

figure(1)
% scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
% hold on
% plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'LineWidth',2)
% 
% plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'LineWidth',2)
% title("Initial and Final Orbits")
% xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
% ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)


options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @multiple_cross_direction_agnostic);
no_event_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

saved_final_state_init_mult_cross = [];
count = 0;

% Multiple Cross
for j = 1:length(xout_lyapunov_init)
% for j = 1:1
    disp(j)
    init_state_0 = [xout_lyapunov_init(j,:), init_mass];
    for i = 1:num_angles
        count = count + 1;

        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction(:,i));
        sol = ode113(fun, [0, 25], init_state_0, options);
        if init_state_0(2) == 0
            seconds_init = sol.xe(3);
            saved_final_state_init_mult_cross(count,:) = [sol.ye(1,3), xout(4,3), xout(5,3), i, j, seconds_init];
        else
            seconds_init = sol.xe(2);
            saved_final_state_init_mult_cross(count,:) = [sol.ye(1,2), xout(4,2), xout(5,2), i, j, seconds_init];
        end

        [tout, xout] = ode113(fun, [0, seconds_init], init_state_0, no_event_opts);
        plot(xout(:,1), xout(:,2), 'Color', 'b')
    end
end

saved_final_state_final_mult_cross = [];
count = 0;

% Multiple Cross
for j = 1:length(xout_lyapunov_final)
    disp(j)
    init_state_0 = [xout_lyapunov_final(j,:), init_mass];
    for i = 1:num_angles
        count = count + 1;

        % Start from pointing straight down and rotate ccw
        fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction(:,i));
        sol = ode113(fun, [0, -25], init_state_0, options);
        if init_state_0(2) == 0
            seconds_final = sol.xe(3);
            saved_final_state_final_mult_cross(count,:) = [sol.ye(1,3), xout(4,3), xout(5,3), i, j, seconds_final];
        else
            seconds_final = sol.xe(2);
            saved_final_state_final_mult_cross(count,:) = [sol.ye(1,2), xout(4,2), xout(5,2), i, j, seconds_final];
        end

        [tout, xout] = ode113(fun, [0, seconds_final], init_state_0, no_event_opts);
        plot(xout(:,1), xout(:,2), 'Color', 'r')
    end
end

hold off
grid on
legend("L2", "Initial Lyapunov Orbit", "Final Lyapunov Orbit")


%%

figure(2)
scatter3(saved_final_state_init_mult_cross(:,1), saved_final_state_init_mult_cross(:,2), saved_final_state_init_mult_cross(:,3), 'filled', 'blue')
% hold on
scatter3(saved_final_state_final_mult_cross(:,1), saved_final_state_final_mult_cross(:,2), saved_final_state_final_mult_cross(:,3), "filled", 'red')
% xlabel("x")
% ylabel("xdot")
% zlabel("ydot")
% legend("Initial Orbit", "Final Orbit")
title("2-Cross Poincar\'e Map", 'Interpreter','latex')
        

%% Plot

figure(3)
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
plot(xout_lyapunov_init(:,1), xout_lyapunov_init(:,2), 'red', 'LineWidth', 2)
plot(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), 'blue', 'LineWidth', 2)
hold off
grid on
title("Initial and Final Orbits")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
legend("L2", "Initial Orbit", "Final Orbit")


%% Plot uncorrected trajectory


uncorrected_init_idx = 1901;
uncorrected_final_idx = 2309;
uncorrected_init_i = saved_final_state_init_mult_cross(uncorrected_init_idx,4);
uncorrected_init_j = saved_final_state_init_mult_cross(uncorrected_init_idx,5);
uncorrected_init_seconds = saved_final_state_init_mult_cross(uncorrected_init_idx,6);
uncorrected_final_i = saved_final_state_final_mult_cross(uncorrected_final_idx,4);
uncorrected_final_j = saved_final_state_final_mult_cross(uncorrected_final_idx,5);
uncorrected_final_seconds = saved_final_state_final_mult_cross(uncorrected_final_idx,6);
uncorrected_init_state0 = xout_lyapunov_init(uncorrected_init_j,:);
uncorrected_final_state0 = xout_lyapunov_final(uncorrected_final_j,:);

fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction(:,uncorrected_init_i));
[tout, xout] = ode113(fun, [0, uncorrected_init_seconds], [uncorrected_init_state0, init_mass]);
x_1_f = xout(end,:)';
figure()
plot(xout(:,1), xout(:,2), 'Color', 'b', 'LineWidth', 2)
hold on

fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction(:,uncorrected_final_i));
[tout, xout] = ode113(fun, [0, uncorrected_final_seconds], [uncorrected_final_state0, init_mass]);
plot(xout(:,1), xout(:,2), 'Color', 'r', 'LineWidth', 2)
hold off
grid on
title("Initial and Final Transfers - Uncorrected")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
legend("Initial Transfer", "Final Transfer")

%% Correction

V0 = [uncorrected_init_state0'; uncorrected_init_seconds; uncorrected_final_state0'; uncorrected_final_seconds];

e1 = x_1_f - uncorrected_init_state0';
F = e1;

V_soln = correction(V0, mu, true);




%% Jacobi constant

function C = jacobi_constant(x, mu)
    C = u_star_times_2(x(1), x(2), x(3), mu) - x(4)^2 - x(5)^2 - x(6)^2;
end

