clear; clc; close all;

% ASEN 5010 - Project, Programming Assignment, Task 7
% Spring 2025
% Jash Bhalavat

%% Given

r_mars = 3396.19; % km
h = 400; % km
r_lmo = r_mars + h;

mu_mars = 42828.3; % km^3/s^2
theta_dot_lmo = sqrt(mu_mars/r_lmo^3);

r_gmo = 20424.2; % km
theta_dot_gmo = 0.0000709003; % rad/s

antenna_n_hat_b = [-1, 0, 0]';
sensor_n_hat_b = [1, 0, 0]';
sa_n_hat_b = [0, 0, 1]';

sigma_BN_t0 = [0.3, -0.4, 0.5]';
omega_BN_B_t0 = [1, 1.75, -2.2]' * pi/180; % rad/s
I_B = diag([10, 5, 7.5]); % kg-m^2

omega_0_lmo = 20 * pi/180;
i_0_lmo = 30 * pi/180;
theta_0_lmo = 60 * pi/180;
period_lmo = 2*pi*sqrt(r_lmo^3/mu_mars);

omega_0_gmo = 0 * pi/180;
i_0_gmo = 0 * pi/180;
theta_0_gmo = 250 * pi/180;
period_gmo = 2*pi*sqrt(r_gmo^3/mu_mars);

n_lmo = ceil(period_lmo); % time steps per orbit
n_gmo = ceil(period_gmo);

t_lmo = 0:n_lmo;
t_gmo = 0:n_gmo;

%% Task 1

for i = 1:length(t_lmo)
    theta_lmo(i,1) = wrapTo2Pi(theta_dot_lmo*t_lmo(i) + theta_0_lmo);
    ea_ON_lmo(i,:) = [omega_0_lmo, i_0_lmo, theta_lmo(i,1)];
    [r_N_lmo(i,:), r_dot_N_lmo(i,:)] = inertial_velocity(r_lmo, ea_ON_lmo(i,:), theta_dot_lmo);
end

for i = 1:length(t_gmo)
    theta_gmo(i,1) = wrapTo2Pi(theta_dot_gmo*t_gmo(i) + theta_0_gmo);
    ea_ON_gmo(i,:) = [omega_0_gmo, i_0_gmo, theta_gmo(i,1)];
    [r_N_gmo(i,:), r_dot_N_gmo(i,:)] = inertial_velocity(r_gmo, ea_ON_gmo(i,:), theta_dot_gmo);
end


%% Task 2

% [NT] = [t1_N, t2_N, t3_N]
% Therefore, [NH] = [ir_N, itheta_N, ih_N];

% @300 seconds, t_lmo index = 301
t_300_ind = 301;

DCM_HN_at_300 = orbit_frame(t_300_ind, r_N_lmo, r_dot_N_lmo);

%% Task 3

DCM_RsN_at_0 = sun_reference_frame();

function DCM_RsN = sun_reference_frame()
    r1_N = [-1, 0, 0]';
    r3_N = [0, 1, 0]';
    r2_N = cross(r3_N, r1_N);
    
    % NR = [r1_N, r2_N, r3_N]
    NRs = [r1_N, r2_N, r3_N];
    DCM_RsN = NRs';
end

%% Task 4

for i = 1:length(t_lmo)
    [DCM_RnN(:,:,i), omega_RnN_N(:,i)] = nadir_reference_frame(r_N_lmo(i,:)', r_dot_N_lmo(i,:)', theta_dot_lmo);
end

DCM_RnN_at_330 = DCM_RnN(:,:,331);
omega_RnN_N_at_330 = omega_RnN_N(:,331);

% DCM_RnN_at_330 = nadir_reference_frame_orientation;

function [DCM_RnN, omega_RnN_N] = nadir_reference_frame(r_N, v_N, theta_dot)
    r1_N = -r_N/norm(r_N);
    r2_N = v_N/norm(v_N);
    r3_N = cross(r1_N, r2_N);
    
    % NR = [r1_N, r2_N, r3_N]
    NRn = [r1_N, r2_N, r3_N];
    DCM_RnN = NRn';

    % omega_Rn/N_Rn = [0, 0, -theta_dot]
    omega_RnN_Rn = [0, 0, -theta_dot]';
    omega_RnN_N = DCM_RnN' * omega_RnN_Rn;
end


%% Task 5

[DCM_RcN_at_330, omega_RcN_N_at_330] = gmo_reference_frame(r_N_lmo(332,:)', r_N_gmo(332,:)', t_lmo(332), r_N_lmo(331,:)', r_N_gmo(331,:)', t_lmo(331));
% [DCM_RcN_at_330, omega_RcN_N_at_330] = gmo_reference_frame(r_N_lmo(16502,:)', r_N_gmo(16502,:)', t_lmo(16502), r_N_lmo(16501,:)', r_N_gmo(16501,:)', t_lmo(16501));

function [DCM_RcN, omega_RcN_N] = gmo_reference_frame(r_LMO_N, r_GMO_N, t, prev_r_LMO_N, prev_r_GMO_N, prev_t)
    % Time step - 330 seconds
    prev_r_GMO_wrt_LMO_N = prev_r_GMO_N - prev_r_LMO_N;
    prev_r1 = -prev_r_GMO_wrt_LMO_N/norm(prev_r_GMO_wrt_LMO_N);
    n3 = [0 0 1]';
    prev_r2 = cross(prev_r_GMO_wrt_LMO_N, n3) / norm(cross(prev_r_GMO_wrt_LMO_N, n3));
    prev_r3 = cross(prev_r1, prev_r2);

    prev_NRc = [prev_r1, prev_r2, prev_r3];
    prev_DCM_RcN = prev_NRc';

    % Time step - 331 seconds
    r_GMO_wrt_LMO_N = r_GMO_N - r_LMO_N;
    r1 = -r_GMO_wrt_LMO_N/norm(r_GMO_wrt_LMO_N);
    n3 = [0 0 1]';
    r2 = cross(r_GMO_wrt_LMO_N, n3) / norm(cross(r_GMO_wrt_LMO_N, n3));
    r3 = cross(r1, r2);

    NRc = [r1, r2, r3];
    DCM_RcN = NRc';

    % To calculate omega_RcN_N, use the following property:
    % C_dot = -omega_tilde * C
    % omega_tilde = -C_dot * C'
    % Where C_dot is numerically calculated by DCM_diff/dt
    dt = t - prev_t;
    dcm_diff = DCM_RcN - prev_DCM_RcN;
    dcm_RcN_dot = dcm_diff/dt;
    omega_tilde_Rc = -dcm_RcN_dot * DCM_RcN';
    omega_RcN_Rc = [-omega_tilde_Rc(2,3); omega_tilde_Rc(1,3); -omega_tilde_Rc(1,2)];
    omega_RcN_N = DCM_RcN' * omega_RcN_Rc;
end



%% Task 6

% At t=0
DCM_RsN_at_0 = sun_reference_frame();
omega_RsN_N_at_0 = [0 0 0]';
[DCM_RnN_at_0, omega_RnN_N_at_0] = nadir_reference_frame(r_N_lmo(1,:)', r_dot_N_lmo(1,:)', theta_dot_lmo);
[DCM_RcN_at_0, omega_RcN_N_at_0] = gmo_reference_frame(r_N_lmo(1,:)', r_N_gmo(1,:)', t_lmo(1), r_N_lmo(2,:)', r_N_gmo(2,:)', t_lmo(2));

[sigma_BRs_at_0, omega_BRs_B_at_0] = attitude_error(t_lmo(1), sigma_BN_t0, omega_BN_B_t0, DCM_RsN_at_0, omega_RsN_N_at_0);
[sigma_BRn_at_0, omega_BRn_B_at_0] = attitude_error(t_lmo(1), sigma_BN_t0, omega_BN_B_t0, DCM_RnN_at_0, omega_RnN_N_at_0);
[sigma_BRc_at_0, omega_BRc_B_at_0] = attitude_error(t_lmo(1), sigma_BN_t0, omega_BN_B_t0, DCM_RcN_at_0, omega_RcN_N_at_0);

function [sigma_BR, omega_BR_B] = attitude_error(t, sigma_BN, omega_BN_B, DCM_RN, omega_RN_N)
    DCM_BN = mrp_to_dcm(sigma_BN);
    DCM_BR = DCM_BN * DCM_RN';
    sigma_BR = dcm_to_mrp(DCM_BR);

    omega_NR_B = DCM_BN * -omega_RN_N;
    omega_BR_B = omega_BN_B + omega_NR_B;
end


%% Task 7

state_0 = [sigma_BN_t0; omega_BN_B_t0];
t = 0:n_lmo;
u_in = [0; 0; 0];

[state, u] = rk4(@(state, t, u_in, I_B)diff_eq(state, t, u_in, I_B), t, state_0, I_B, u_in);

% t = 500 s @ t ind --> 1001
omega_BN_B_at_500 = state(501, 4:6)';
H_B_at_500 = I_B * omega_BN_B_at_500;
print_array(H_B_at_500)

T_at_500 = 1/2 * omega_BN_B_at_500' * I_B * omega_BN_B_at_500;
print_array(T_at_500)

sigma_BN_at_500 = state(501, 1:3)';
print_array(sigma_BN_at_500)
DCM_BN_at_500 = mrp_to_dcm(sigma_BN_at_500);
omega_BN_N_at_500 = DCM_BN_at_500' * omega_BN_B_at_500;
I_N = DCM_BN_at_500' * I_B * DCM_BN_at_500;
H_N_at_500 = I_N * omega_BN_N_at_500;
print_array(H_N_at_500)

u_in = [0.01; -0.01; 0.02];
[state_control, u] = rk4(@(state, t, u_in, I_B)diff_eq(state, t, u_in, I_B), t, state_0, I_B, u_in);

sigma_BN_at_100 = state_control(101, 1:3)';
print_array(sigma_BN_at_100)


for i = 1:length(t)
    omega_BN_B = state(i, 4:6)';
    H_B(:,i) =  I_B * omega_BN_B;
    T(i) = 1/2 * omega_BN_B' * I_B * omega_BN_B;
    sigma_BN(:,i) = state(i, 1:3)';
    DCM_BN = mrp_to_dcm(sigma_BN(:,i));
    omega_BN_N(:,i) = DCM_BN' * omega_BN_B;
    I_N = DCM_BN' * I_B * DCM_BN;
    H_N(:,i) = I_N * omega_BN_N(:,i);
end


%% Plot

figure()
subplot(5, 1, 1)
plot(t, H_B, 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("$^BH [Nms]$", Interpreter="latex")
legend("$^BH(1)$", "$^BH(2)$", "$^BH(3)$", Interpreter='latex')
grid on
title("Angular Momentum in Body Frame")
xlim([0, 500])

subplot(5, 1, 2)
plot(t, T, 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("$T [J]$", Interpreter="latex")
grid on
title("Rotational Kinetic Energy")
xlim([0, 500])

subplot(5, 1, 3)
plot(t, sigma_BN, 'LineWidth', 2)
legend("$\sigma_{B/N}(1)$", "$\sigma_{B/N}(2)$", "$\sigma_{B/N}(3)$", 'FontSize', 12, Interpreter="latex")
ylabel("$\sigma_{B/N}$", 'FontSize', 14, Interpreter="latex")
xlabel("Time [sec]")
grid on
title("Inertial to Body Frame MRP")
xlim([0, 500])

subplot(5, 1, 4)
plot(t, omega_BN_N, 'LineWidth',2)
legend("$^N\omega_{B/N}(1)$", "$^N\omega_{B/N}(2)$", "$^N\omega_{B/N}(3)$", 'FontSize', 12, Interpreter="latex")
ylabel("$^N\omega_{B/N}$ [rad/sec]", 'FontSize', 14, Interpreter="latex")
xlabel("Time [sec]")
grid on
title("Angular Velocity of Body w.r.t. Inertial Frame in Inertial Frame")
xlim([0, 500])

subplot(5, 1, 5)
plot(t, H_N, 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("$^NH [Nms]$", Interpreter="latex")
legend("$^NH(1)$", "$^NH(2)$", "$^NH(3)$", Interpreter='latex')
grid on
title("Angular Momentum in Inertial Frame")
xlim([0, 500])

sgtitle("Task 7 - Outputs")


%% Functions

function state_dot = diff_eq(state, t, u, I)
    % 6x1 state EOMs
    % Inputs:
    % state - [sigma_BN'; omega_BN_B']
    % I - 3x3 principal MOI [kg-m^2]
    % 
    % Output:
    % state_dot - [sigma_BN_dot', omega_BN_B_dot']
    sigma_BN = state(1:3)';
    omega_BN_B = state(4:6)';

    sigma_sq = dot(sigma_BN, sigma_BN);
    sigma_BN_dot = 1/4 * ((1-sigma_sq)*eye(3) + 2*skew_symmetric(sigma_BN) + 2*(sigma_BN*sigma_BN')) * omega_BN_B;

    omega_BN_B_dot(1,1) = -(I(3,3) - I(2,2))/I(1,1) * omega_BN_B(2) * omega_BN_B(3) + u(1)/I(1,1);
    omega_BN_B_dot(2,1) = -(I(1,1) - I(3,3))/I(2,2) * omega_BN_B(1) * omega_BN_B(3) + u(2)/I(2,2);
    omega_BN_B_dot(3,1) = -(I(2,2) - I(1,1))/I(3,3) * omega_BN_B(2) * omega_BN_B(1) + u(3)/I(3,3);
    
    state_dot = [sigma_BN_dot; omega_BN_B_dot];
end


function dcm = mrp_to_dcm(mrp)
    mrp_sq = norm(mrp)^2;
    mrp_tilde = skew_symmetric(mrp);
    dcm = eye(3) + (8 * mrp_tilde * mrp_tilde - 4 * (1-mrp_sq) * mrp_tilde)/(1+mrp_sq)^2;
end

function mrp = dcm_to_mrp(dcm)
    zeta = sqrt(trace(dcm) + 1);
    mrp = 1/(zeta*(zeta+2)) * [dcm(2,3)-dcm(3,2); dcm(3,1)-dcm(1,3); dcm(1,2)-dcm(2,1)];
end

function DCM_HN = orbit_frame(t, r_vec, r_dot_vec)
    r = r_vec(t,:)';
    r_dot = r_dot_vec(t,:)';
    i_r_N = r / norm(r);
    i_h_N = cross(r, r_dot) / norm(cross(r, r_dot));
    i_theta_N = cross(i_h_N, i_r_N);
    NH = [i_r_N, i_theta_N, i_h_N];
    DCM_HN = NH';
end


function [r_N, r_dot_N] = inertial_velocity(r, ea, theta_dot)
    % Inertial s/c velocity vector for a circular orbit
    % Inputs
    % r - radius
    % ea - 3-1-3 euler angles
    % 
    % Outputs
    % r_N - inertial position
    % r_dot_N - inertial velocity

    r_O = [r, 0, 0]';
    omega_ON_O = [0, 0, theta_dot]';
    
    ON = R3(ea(3)) * R1(ea(2)) * R3(ea(1));


    r_N = ON' * r_O;
    r_dot_N = ON' * cross(omega_ON_O, r_O);

end




