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

lyapunov_orbits = load("V_family_L1_Lyapunov.mat").V_family;
halo_orbits = load("V_family_halo_L1_northern.mat").V_family;
em_eq_pts = load("eq_points.mat").em_eq_pts;

l2_pos = [em_eq_pts(2,:), 0];
l1_pos = [em_eq_pts(1,:), 0];
p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure(1)
scatter3(l1_pos(1), l1_pos(2), 0, 'filled', 'red')
hold on

% Set options for ode113
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[tout_lyapunov_init, xout_lyapunov_final] = ode113(@(t,state)CR3BP(state, mu), [0, lyapunov_orbits(7,5)], lyapunov_orbits(1:6,5), options);
plot3(xout_lyapunov_final(:,1), xout_lyapunov_final(:,2), xout_lyapunov_final(:,3), 'LineWidth',2)

[tout_halo_final, xout_halo_init] = ode113(@(t,state)CR3BP(state, mu), [0, halo_orbits(7,200)], halo_orbits(1:6,200), options);
plot3(xout_halo_init(:,1), xout_halo_init(:,2), xout_halo_init(:,3), 'LineWidth',2)

title("Initial and Final Orbits")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)

%% Transfer

% Initial orbit, initial state
init_mass = 10;

function [value, isterminal, direction] = zero_y_hyperplane_init(t, state)
    value = state(2);          % Stop when y = 0
    isterminal = 1;     % 1 = stop integration
    direction = 1;     % Direction agnostic
end

function [value, isterminal, direction] = l1_x_hyperplane_init(t, state)
    value = state(1)-8.369151317503717e-01; % Stop when x is @ l1
    isterminal = 1;     % 1 = stop integration
    direction = 0;     % Direction agnostic
end

function [value, isterminal, direction] = l1_x_hyperplane_final(t, state)
    value = state(1)-8.369151317503717e-01; % Stop when x is @ l1
    isterminal = 1;     % 1 = stop integration
    direction = 0;     % Direction agnostic
end

function [value, isterminal, direction] = zero_y_hyperplane_final(t, state)
    value = state(2);          % Stop when y = 0
    isterminal = 1;     % 1 = stop integration
    direction = 0;     % Direction agnostic
end

num_angles = 5;
angle1 = linspace(0, pi, num_angles);
angle2 = linspace(0, pi, num_angles);
angle3 = linspace(0, pi, num_angles);

% Straight right vector
neg_y_vec = [1, 0, 0]';

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @l1_x_hyperplane_init);

saved_final_state_init = [];
count = 0;

% for j = 1:2:length(xout_halo_init)
%     init_state_0 = [xout_halo_init(j,:), init_mass];
%     disp(j)
%     for i = 1:num_angles
%         for k = 1:num_angles
%             for l = 1:num_angles
%                 dcm = R1(angle1(i))*R2(angle2(k))*R3(angle3(l));
%                 thrust_direction = dcm * neg_y_vec;
%                 fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction);
%                 [tout, xout] = ode113(fun, [0, 10], init_state_0, options);
% 
%                 % Only plot if last point is on the hyperplane beyond the L2 point
%                 count = count + 1;
%                 saved_final_state_init(count,:) = [xout(end,1), xout(end,4), xout(end,5)];
%                 plot3(xout(:,1), xout(:,2), xout(:,3), 'Color', 'b')
%             end
%         end
%     end
% end

for j = 1:2:length(xout_halo_init)
    init_state_0 = [xout_halo_init(j,:), init_mass];
    % Straight down thrust direction
    thrust_direction = [0, 0, 1];
    fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction);
    [tout, xout] = ode113(fun, [0, 10], init_state_0, options);

    % Only plot if last point is on the hyperplane beyond the L2 point
    if abs(xout(end,1)-l1_pos(1)) < 1e-6
        count = count + 1;
        saved_final_state_init(count,:) = xout(end,:);
        plot3(xout(:,1), xout(:,2), xout(:,3), 'Color', 'b')
    end
end

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @l1_x_hyperplane_final);

saved_final_state_final = [];
count = 0;

% Going back from final orbit
% for j = 1:2:length(xout_lyapunov_final)
%     init_state_0 = [xout_lyapunov_final(j,:), init_mass];
%     for i = 1:num_angles
%         % Start from pointing straight down and rotate ccw
%         dcm = R3(angles(i));
%         thrust_direction = dcm * neg_y_vec;
%         fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction);
%         [tout, xout] = ode113(fun, [0, -100], init_state_0, options);
% 
%         % plot(xout(:,1), xout(:,2), xout(:,3), 'Color', 'b')
%         % Only plot is last point is on the hyperplane beyond the L2 point
%         if xout(end,1) > l2_pos(1)
%             count = count + 1;
%             saved_final_state_final(count,:) = [xout(end,1), xout(end,4), xout(end,5)];
%             plot(xout(:,1), xout(:,2), 'Color', 'r')
%         end
%     end
% end

for j = 1:2:length(xout_lyapunov_final)
    init_state_0 = [xout_lyapunov_final(j,:), init_mass];
    thrust_direction = [0, 0, 1];
    fun = @(t,state)CR3BP_with_mass(state, mu, T, Isp, T_star, thrust_direction);
    [tout, xout] = ode113(fun, [0, -100], init_state_0, options);

    % plot(xout(:,1), xout(:,2), xout(:,3), 'Color', 'b')
    % Only plot is last point is on the hyperplane beyond the L2 point
    if abs(xout(end,1)-l1_pos(1)) < 1e-6
        count = count + 1;
        saved_final_state_final(count,:) = xout(end,:);
        plot3(xout(:,1), xout(:,2), xout(:,3), 'Color', 'r')
    end
end

hold off
grid on
legend("L2", "Initial Lyapunov Orbit", "Final Lyapunov Orbit")

%%

figure(2)
% scatter3(saved_final_state_init(:,1), saved_final_state_init(:,2), saved_final_state_init(:,3), 'filled', 'blue')
quiver3(saved_final_state_init(:,1), saved_final_state_init(:,2), saved_final_state_init(:,3), saved_final_state_init(:,4), saved_final_state_init(:,5), saved_final_state_init(:,6), 'blue', 'LineWidth', 2, 'MaxHeadSize', 2)
hold on
% scatter3(saved_final_state_final(:,1), saved_final_state_final(:,2), saved_final_state_final(:,3), "filled", 'red')
quiver3(saved_final_state_final(:,1), saved_final_state_final(:,2), saved_final_state_final(:,3), saved_final_state_final(:,4), saved_final_state_final(:,5), saved_final_state_final(:,6), 'red', 'LineWidth', 2, 'MaxHeadSize', 2)
xlabel("x")
ylabel("y")
zlabel("z")
title("Hyperplane Crossings")






