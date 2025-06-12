clear; clc; close all;

% ASEN 5010 - HW 3, Problem 4
% Spring 2025
% Jash Bhalavat

mass = 1; % kg
v_CN_P_0 = [0, 0, 0]'; % m/s
v_CN_N_0 = v_CN_P_0;
omega_PN_P_0 = [0, 0, 0.85]'; % rad/sec
euler_PN_0 = [0, 0, 0]'; % rad
g = 9.81; % km/s2
I_c_P = [32.5, 32.5, 5]'; % kg-m2
r_NozzleC_P = [0, 0, -3]'; % m
r_NoseC_P = [0, 0, 4]'; % m

% Guessing initial position of rocket CoM wrt N frame
% Letting z be 3 assuming nozzle touches the ground
r_CN_N_0 = [0, 0, 3]'; % [m]

tspan = [0, 60]; % sec

% Thrust at rear nozzle
% thrust_cm = cos(thruster_misalignment) * thruster_mag; % [N] Passes CoM
% thrust_perp = sin(thruster_misalignment) * thruster_mag; % [N]

function thrust_P = thrust(time)
    thruster_misalignment = deg2rad(2.5); % rad
    thruster_mag = 15; % N
    thruster_duration = 10; % sec
    thrust_parallel = 0; % [N] Passes CoM
    thrust_perp = 0; % [N]
    if time <= thruster_duration
        thrust_parallel = thruster_mag * cos(thruster_misalignment);
        thrust_perp = thruster_mag * sin(thruster_misalignment);
    end
    thrust_P = [0, -thrust_perp, thrust_parallel]';
end

%% Part a

[tout_omega_PN_P, omega_PN_P_out] = ode45(@(t, omega)rot_eom(t, omega, I_c_P, r_NozzleC_P), [tspan(1), tspan(2)], omega_PN_P_0);

figure()
subplot(3,1,1)
plot(tout_omega_PN_P, omega_PN_P_out(:,1), 'LineWidth',2)
ylabel("omega_{PN}^P(1) [rad/s]")
grid on

subplot(3,1,2)
plot(tout_omega_PN_P, omega_PN_P_out(:,2), 'LineWidth',2)
ylabel("omega_{PN}^P(2) [rad/s]")
grid on

subplot(3,1,3)
plot(tout_omega_PN_P, omega_PN_P_out(:,3), 'LineWidth',2)
ylabel("omega_{PN}^P(3) [rad/s]")
xlabel("Time [sec]")
sgtitle("Angular velocity of rocket in principal frame")
grid on

function out = L_c_P(time, r)
    % Torque at CoM in P frame
    thrust_at_time_P = thrust(time);
    % Already negative, so don't need to add (-)
    % out = [3*thrust_at_time_P(2), 0, 0]';
    out = cross(r, thrust_at_time_P);
end

function omega_PN_P_dot = rot_eom(time, omega_PN_P, I, r_NozzleC_P)
    L_c_P_at_time = L_c_P(time, r_NozzleC_P);
    omega_PN_P_dot(1,1) = -1/I(1) * (I(3) - I(2)) * omega_PN_P(2)*omega_PN_P(3) + L_c_P_at_time(1);
    omega_PN_P_dot(2,1) = -1/I(2) * (I(1) - I(3)) * omega_PN_P(1)*omega_PN_P(3) + L_c_P_at_time(2);
    omega_PN_P_dot(3,1) = -1/I(3) * (I(2) - I(1)) * omega_PN_P(1)*omega_PN_P(2) + L_c_P_at_time(3);
end

%% Part b

[tout_v_CN_P, v_CN_P_out] = ode45(@(t, v)trans_eom(t, v, mass, g), [tspan(1), tspan(2)], v_CN_P_0);

figure()
subplot(3,1,1)
plot(tout_v_CN_P, v_CN_P_out(:,1), 'LineWidth',2)
ylabel("v_{CN}^P(1) [m/s]")
grid on

subplot(3,1,2)
plot(tout_v_CN_P, v_CN_P_out(:,2), 'LineWidth',2)
ylabel("v_{CN}^P(2) [m/s]")
grid on

subplot(3,1,3)
plot(tout_v_CN_P, v_CN_P_out(:,3), 'LineWidth',2)
ylabel("v_{CN}^P(3) [m/s]")
xlabel("Time [sec]")
grid on
sgtitle("Translational velocity of rocket's CoM in principal frame")

function v_CN_P_dot = trans_eom(time, v_CN_P, mass, g)
    thrust_at_time_P = thrust(time);
    v_CN_P_dot(1,1) = 0;
    v_CN_P_dot(2,1) = thrust_at_time_P(2)/mass;
    v_CN_P_dot(3,1) = (thrust_at_time_P(3) - mass*g)/mass;
end

%% Part c

omega_PN_P_euler_PN_0 = [omega_PN_P_0; euler_PN_0];

[tout_omega_PN_P_euler_PN, omega_PN_P_euler_PN_out] = ode45(@(t, rot_euler)rot_euler_eom(t, rot_euler, I_c_P, r_NozzleC_P), [tspan(1), tspan(2)], omega_PN_P_euler_PN_0);

figure()
subplot(2,1,1)
plot(tout_omega_PN_P_euler_PN, omega_PN_P_euler_PN_out(:,4), 'LineWidth',2)
ylabel("\alpha [rad]")
grid on

subplot(2,1,2)
plot(tout_omega_PN_P_euler_PN, omega_PN_P_euler_PN_out(:,5), 'LineWidth',2)
xlabel("Time [sec]")
ylabel("\beta [rad]")
sgtitle("\alpha and \beta Euler angles")
grid on

function rot_euler_dot = rot_euler_eom(time, rot_euler, I, r_NozzleC_P)
    % rot - 3x1 omega_PN_P
    % euler - 3x1 euler_PN
    omega_PN_P = rot_euler(1:3);
    euler_PN = rot_euler(4:6);
    L_c_P_at_time = L_c_P(time, r_NozzleC_P);
    omega_PN_P_dot(1,1) = -1/I(1) * (I(3) - I(2)) * omega_PN_P(2)*omega_PN_P(3) + L_c_P_at_time(1);
    omega_PN_P_dot(2,1) = -1/I(2) * (I(1) - I(3)) * omega_PN_P(1)*omega_PN_P(3) + L_c_P_at_time(2);
    omega_PN_P_dot(3,1) = -1/I(3) * (I(2) - I(1)) * omega_PN_P(1)*omega_PN_P(2) + L_c_P_at_time(3);

    B_theta = 1/cos(euler_PN(2)) * [cos(euler_PN(3)), -sin(euler_PN(3)), 0;
                cos(euler_PN(2))*sin(euler_PN(3)), cos(euler_PN(2))*cos(euler_PN(3)), 0;
                -sin(euler_PN(2))*cos(euler_PN(3)), sin(euler_PN(2))*sin(euler_PN(3)), cos(euler_PN(2))];
    euler_PN_dot = B_theta * omega_PN_P;
    rot_euler_dot = [omega_PN_P_dot; euler_PN_dot];
end

%% Part d

state0 = [r_CN_N_0; v_CN_N_0; omega_PN_P_0; euler_PN_0];

[tout_state, state_out] = ode45(@(t, state)state_eom(t, state, mass, g, I_c_P, r_NozzleC_P), [tspan(1), tspan(2)], state0);

figure()
plot3(state_out(:,1), state_out(:,2), state_out(:,3), 'LineWidth',2)
xlabel('$$\hat{n}_{1}$$ [m]','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{n}_{2}$$ [m]','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{n}_{3}$$ [m]','Interpreter','Latex', 'FontSize',18)
hold on
scatter3(state_out(1,1), state_out(1,2), state_out(1,3), 'filled', 'black')
scatter3(state_out(end,1), state_out(end,2), state_out(end,3), 'filled', 'blue')
legend("CoM", "CoM @ t=0 sec", "CoM @ t=60 sec")
grid on
title("Center of Mass in 3D inertial space")

figure()
subplot(3,1,1)
plot(tout_state, state_out(:,1), 'LineWidth',2)
ylabel("r_{CN}^N(1) [m]")
grid on
subplot(3,1,2)
plot(tout_state, state_out(:,2), 'LineWidth',2)
ylabel("r_{CN}^N(2) [m]")
grid on
subplot(3,1,3)
plot(tout_state, state_out(:,3), 'LineWidth',2)
ylabel("r_{CN}^N(3) [m]")
xlabel("Time [sec]")
grid on
sgtitle("Center of Mass in inertial space")

function state_dot = state_eom(time, state, mass, g, I, r_NozzleC_P)
    % state - [x_N, y_N, z_N, x_N_dot, y_N_dot, z_N_dot, omega_PN_P(1), omega_PN_P(2),
    % omega_PN_P(3), alpha, beta, gamma]'
    % state_dot - [x_N_dot, y_N_dot, z_N_dot, x_N_dotdot, y_N_dotdot, z_N_dotdot, omega_PN_P(1), omega_PN_P(2),
    % omega_PN_P(3), alpha, beta, gamma]'

    r_CN_N = state(1:3);
    v_CN_N = state(4:6);
    omega_PN_P = state(7:9);
    euler_PN = state(10:12);

    PN = R3(euler_PN(3))*R2(euler_PN(2))*R1(euler_PN(1));

    thrust_at_time_P = thrust(time);
    v_CN_P_dot(1,1) = 0;
    v_CN_P_dot(2,1) = thrust_at_time_P(2)/mass;
    v_CN_P_dot(3,1) = (thrust_at_time_P(3) - mass*g)/mass;

    r_CN_N_dot = v_CN_N;
    v_CN_N_dot = PN' * v_CN_P_dot;

    L_c_P_at_time = L_c_P(time, r_NozzleC_P);
    omega_PN_P_dot(1,1) = -1/I(1) * (I(3) - I(2)) * omega_PN_P(2)*omega_PN_P(3) + L_c_P_at_time(1);
    omega_PN_P_dot(2,1) = -1/I(2) * (I(1) - I(3)) * omega_PN_P(1)*omega_PN_P(3) + L_c_P_at_time(2);
    omega_PN_P_dot(3,1) = -1/I(3) * (I(2) - I(1)) * omega_PN_P(1)*omega_PN_P(2) + L_c_P_at_time(3);

    B_theta = 1/cos(euler_PN(2)) * [cos(euler_PN(3)), -sin(euler_PN(3)), 0;
                cos(euler_PN(2))*sin(euler_PN(3)), cos(euler_PN(2))*cos(euler_PN(3)), 0;
                -sin(euler_PN(2))*cos(euler_PN(3)), sin(euler_PN(2))*sin(euler_PN(3)), cos(euler_PN(2))];
    euler_PN_dot = B_theta * omega_PN_P;

    state_dot = [r_CN_N_dot; v_CN_N_dot; omega_PN_P_dot; euler_PN_dot];
end


%% Part e

euler_PN = state_out(:,10:12);
r_CN_N = state_out(:,1:3);

for i = 1:length(tout_state)
    PN = R3(euler_PN(i, 3))*R2(euler_PN(i, 2))*R1(euler_PN(i, 1));
    r_NoseN_N(i,:) = (r_CN_N(i,:)' + PN' * r_NoseC_P)';
end

figure()
plot(r_NoseN_N(:,1), r_NoseN_N(:,2), 'LineWidth',2)
hold on
scatter(r_NoseN_N(1,1), r_NoseN_N(1,2), 'filled', 'black')
scatter(r_NoseN_N(end,1), r_NoseN_N(end,2), 'filled', 'blue')
grid on
legend("Rocket Nose", "Starting Point", "Final Point")
xlabel('$$\hat{n}_{1}$$ [m]','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{n}_{2}$$ [m]','Interpreter','Latex', 'FontSize',18)
title("Projection of rocket's nose in inertial XY plane")

%% Extra

% function euler_PN_dot = euler_eom(time, euler_PN, omega_PN_P)
%     B_theta = 1/cos(euler_PN(2)) * [cos(euler_PN(3)), -sin(euler_PN(3)), 0;
%                 cos(euler_PN(2))*sin(euler_PN(3)), cos(euler_PN(2))*cos(euler_PN(3)), 0;
%                 -sin(euler_PN(2))*cos(euler_PN(3)), sin(euler_PN(2))*sin(euler_PN(3)), cos(euler_PN(2))];
%     euler_PN_dot = B_theta * omega_PN_P;
% end



