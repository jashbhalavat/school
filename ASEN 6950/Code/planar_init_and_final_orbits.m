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

num_angles = 50;
angles = linspace(0, 2*pi, num_angles);

% Straight down vector
neg_y_vec = [0, -1, 0]';

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @zero_y_hyperplane_init, 'OutputFcn', @countNegCrossings);

saved_final_state_init = [];
count = 0;


for j = 1:length(xout_lyapunov_init)
% for j = 1:2
    disp(j)
    init_state_0 = [xout_lyapunov_init(j,:), init_mass];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        dcm = R3(angles(i));
        thrust_direction = dcm * neg_y_vec;
        fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction);
        % clear zero_y_hyperplane_init;
        clear countNegCrossings;
        [tout, xout] = ode113(fun, [0, 100], init_state_0, options);

        % plot(xout(:,1), xout(:,2), xout(:,3), 'Color', 'b')
        % Only plot if last point is on the hyperplane beyond the L2 point
        % if xout(end,1) > l2_pos(1)
        if xout(end,1) >= lyapunov_final_x_min && xout(end,1) <= lyapunov_final_x_max
            count = count + 1;
            saved_final_state_init(count,:) = [xout(end,1), xout(end,4), xout(end,5)];
            plot(xout(:,1), xout(:,2), 'Color', 'b')
        end
    end
end

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @zero_y_hyperplane_final, 'OutputFcn', @countPosCrossings);

saved_final_state_final = [];
count = 0;

% Going back from final orbit
for j = 1:length(xout_lyapunov_final)
% for j = 1:2
    disp(j)
    init_state_0 = [xout_lyapunov_final(j,:), init_mass];
    for i = 1:num_angles
        % Start from pointing straight down and rotate ccw
        dcm = R3(angles(i));
        thrust_direction = dcm * neg_y_vec;
        fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction);
        % clear zero_y_hyperplane_final;
        clear countPosCrossings;
        [tout, xout] = ode113(fun, [0, -100], init_state_0, options);

        % plot(xout(:,1), xout(:,2), xout(:,3), 'Color', 'b')
        % Only plot is last point is on the hyperplane beyond the L2 point
        % if xout(end,1) > l2_pos(1)
        if xout(end,1) >= lyapunov_final_x_min && xout(end,1) <= lyapunov_final_x_max
            count = count + 1;
            saved_final_state_final(count,:) = [xout(end,1), xout(end,4), xout(end,5)];
            plot(xout(:,1), xout(:,2), 'Color', 'r')
        end
    end
end

hold off
grid on
legend("L2", "Initial Lyapunov Orbit", "Final Lyapunov Orbit")

%%

figure(2)
scatter3(saved_final_state_init(:,1), saved_final_state_init(:,2), saved_final_state_init(:,3), 'filled', 'blue')
hold on
scatter3(saved_final_state_final(:,1), saved_final_state_final(:,2), saved_final_state_final(:,3), "filled", 'red')
xlabel("x")
ylabel("xdot")
zlabel("ydot")
legend("Initial Orbit", "Final Orbit")

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



%% Jacobi constant

function C = jacobi_constant(x, mu)
    C = u_star_times_2(x(1), x(2), x(3), mu) - x(4)^2 - x(5)^2 - x(6)^2;
end

