clear; clc; close all;

% ASEN 5044 - Final Project - Part 1 - Non-Linear 
% Jash Bhalavat
% Fall 2024

% Constants
mu = 398600; % km^3/s^2
r0 = 6678; % km
n = 4; % state variables
m = 2; % control variables
R_E = 6378; % km, earth radius
omega_E = (2*pi)/86400; % rad/s
h = 300; % km, satellite altitude
dt = 10; % s



%% ODE45 NL Dynamics
xnom_0 = [r0; 0; 0; r0*sqrt(mu/r0^3)];
dx0 = [0; 0.075; 0; -0.021];

t = 0:dt:14000;

xdot = @(x, u) [x(2);
    -mu*x(1)/((x(1)^2 + x(3)^2)^(3/2));
    x(4);
    -mu*x(3)/((x(1)^2 + x(3)^2)^(3/2))];

x0 = xnom_0 + dx0;

xdotwrap = @(t, x) xdot(x, [0;0]);

options = odeset('RelTol',1e-12);
[tout, xnl] = ode45(xdotwrap, t, x0, options);
xnl = xnl';


figure();
subplot(4, 1, 1);
plot(t, xnl(1, :));
title('Total States vs. Time, ODE45');
ylabel('X (km)');
grid on;
subplot(4, 1, 2);
plot(t, xnl(2, :));
ylabel('X-dot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, xnl(3, :));
ylabel('Y (km)');
grid on;
subplot(4, 1, 4);
plot(t, xnl(4, :));
ylabel('Y-dot (km/s)');
xlabel('Time t (s)');
grid on;

%% ODE45 NL Measurements

% Ground station initial location
number_of_ground_stations = 12;
for i = 1:number_of_ground_stations
    theta_i_0(i, 1) = (i-1)*pi/6;
    X_i_0(i, 1) = R_E * cos(theta_i_0(i));
    Y_i_0(i, 1) = R_E * sin(theta_i_0(i));
    X_i_0_dot(i, 1) = - R_E * omega_E * sin(theta_i_0(i));
    Y_i_0_dot(i, 1) = R_E * omega_E * cos(theta_i_0(i));
end

% Store in array
theta_i = [theta_i_0];
X_i = [X_i_0];
Y_i = [Y_i_0];
X_i_dot = [X_i_0_dot];
Y_i_dot = [Y_i_0_dot];

% The way atan2/atan calculate the angles are not matching with the
% problem (for ex - 91° is calculated to be -89° and while that can be
% dealt with some logic statements, it's way easier to just add the
% first value calculated because the 1st GS starts at 0 and all the
% ground stations increase the same amount at each time step.
theta_i_increase = omega_E * dt;

% Calculate ground station positions for each time step
% Manually taking derivative of X_i and Y_i
% X_i = R_E * cos(omega_E*t + theta_i_0)
% So, X_i_dot = - R_E * sin(omega_E*t + theta_i_0) * omega_E
% Y_i = R_E * sin(omega_E*t + theta_i_0)
% So, Y_i_dot = R_E * cos(omega_E*t + theta_i_0) * omega_E
for time = 2:length(t)
    for station = 1:number_of_ground_stations
        X_i(station, time) = R_E * cos(omega_E*t(time) + theta_i_0(station));
        Y_i(station, time) = R_E * sin(omega_E*t(time) + theta_i_0(station));
        X_i_dot(station, time) = - R_E * omega_E * sin(omega_E*t(time) + theta_i_0(station));
        Y_i_dot(station, time) = R_E * omega_E * cos(omega_E*t(time) + theta_i_0(station));

        % theta_i(station, time) = theta_i(station, time-1) + theta_i_increase;
        % if theta_i(station, time) > 2*pi
        %     theta_i(station, time) = theta_i(station, time) - 2*pi;
        % end
        theta_i(station, time) = atan2(Y_i(time), X_i(time));
    end
end

% GS visibility elevation range
ground_station_elevation_range_lower = theta_i - pi/2;
ground_station_elevation_range_upper = theta_i + pi/2;

% This array keeps track of which GS is generating data at all points. If a
% GS is within FOV, that station number will be stored (1-12) otherwise
% it's set to NaN
p_i = NaN(number_of_ground_stations, length(t));

% Go through time and find range, range rate, and elevation using nonlinear
% dynamics states
for time = 1:length(t)
    x_t = xnl(1, time);
    x_dot_t = xnl(2, time);
    y_t = xnl(3, time);
    y_dot_t = xnl(4, time);
    
    counter = 1;
    
    for station = 1:number_of_ground_stations
        rho_i(station, time) = sqrt((x_t - X_i(station, time))^2 + (y_t - Y_i(station, time))^2);
        rho_dot_i(station, time) = ((x_t - X_i(station, time))*(x_dot_t - X_i_dot(station, time)) + (y_t - Y_i(station, time))*(y_dot_t - Y_i_dot(station, time))) / rho_i(station, time);
        
        phi_i(station, time) = atan2(y_t - Y_i(station, time), x_t - X_i(station, time));

        sc_from_gs = [x_t - X_i(station, time), y_t - Y_i(station, time)];
        angle = ang_diff(sc_from_gs, [X_i(station, time), Y_i(station, time)]);
        if abs(angle) <= pi/2
            p_i(counter, time) = station;
            counter = counter + 1;
        end


        % if phi_i(station, time) >= ground_station_elevation_range_lower(station, time) && phi_i(station, time) <= ground_station_elevation_range_upper(station, time)
        %     p_i(counter, time) = station;
        %     counter = counter + 1;
        % end
    end
end
%% 
count_12 = ones(12, 1);
time_id_ground_station = zeros(12, length(t));
for i = 1:length(t)
    col = p_i(:, i);
    for j = 1:number_of_ground_stations
        if ~isnan(col(j))
            time_id_ground_station(col(j), count_12(col(j))) = i;
            count_12(col(j)) = count_12(col(j)) + 1;
        end
    end
end

colors_12 = {'red', 'green', 'blue', 'cyan', 'yellow', 'black', '#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#A2142F'};
figure()
for i = 1:number_of_ground_stations
    indices_gs = time_id_ground_station(i, 1:count_12(i)-1);
    scatter(t(indices_gs), rho_i(i, indices_gs), 'LineWidth', 1.25, 'Marker', 'x', 'MarkerEdgeColor', colors_12{i})
    hold on
end
hold off

figure()
for i = 1:number_of_ground_stations
    indices_gs = time_id_ground_station(i, 1:count_12(i)-1);
    scatter(t(indices_gs), rho_dot_i(i, indices_gs), 'LineWidth', 1.25, 'MarkerEdgeColor', colors_12{i})
    hold on
end
hold off

figure()
for i = 1:number_of_ground_stations
    indices_gs = time_id_ground_station(i, 1:count_12(i)-1);
    scatter(t(indices_gs), phi_i(i, indices_gs), 'LineWidth', 1.5, 'MarkerEdgeColor', colors_12{i})
    hold on
end
hold off


function angle = ang_diff(V1, V2)
    angle = acos(dot(V1, V2) / (norm(V1) * norm(V2)));
end



