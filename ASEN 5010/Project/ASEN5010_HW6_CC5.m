clear; clc; close all;

% ASEN 5010 - HW 6, CC1
% Spring 2025
% Jash Bhalavat

I_B = diag([100, 75, 80]); % kg-m2
K = 5;
P_i = sqrt(K*I_B);
% gainList = [22.360679774997898, 19.364916731037084, 20 ] # P1, P2, P3

T_i = 2*I_B/P_i;
% decayTimes = [8.944271909999159, 7.745966692414834, 8.0 ] # seconds, T1, T2, T3

sigma_BN_t0 = [0.1, 0.2, -0.1]';
omega_BN_B_t0 = [30, 10, -20]' * pi/180;

% P = 10 * eye(3);
t_span = 0:120;

state_0 = [sigma_BN_t0; omega_BN_B_t0];

[state, u] = rk4(@(state, t, u_B, I_B)diff_eq(state, t, u_B, I_B), t_span, state_0, I_B, K, P_i);

sigma_BN_at_30 = state(31, 1:3);
sigma_BN_norm_at_30 = sqrt(sigma_BN_at_30(1)^2 + sigma_BN_at_30(2)^2 + sigma_BN_at_30(3)^2);

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

function [x, u_B] = rk4(fn, t, x0, I_B, K, P)

    dt = t(2) - t(1);
    
    % x is nx6
    x = x0';

    for i = 1:length(t)
        sigma_BN = x(i,1:3)';
        omega_BN_B = x(i,4:6)';
        % [DCM_RsN, omega_RsN_N] = sun_reference_frame();
        % [sigma_BRs(:,i), omega_BRs_B(:,i)] = attitude_error(t(i), sigma_BN, omega_BN_B, DCM_RsN, omega_RsN_N);

        % sigma_RN(:,i) = zeros([3,1]);
        omega_RN_B = zeros([3,1]);
        L = zeros([3,1]);
        omega_dot_RN_B = zeros([3,1]);
        
        u_B(i,:) = (-K*sigma_BN -P*omega_BN_B + I_B*(omega_dot_RN_B - cross(omega_BN_B, omega_RN_B)) + skew_symmetric(omega_BN_B)*I_B*omega_BN_B - L)';

        k1 = dt * fn(x(i,:), t(i), u_B(i,:), I_B)';
        k2 = dt * fn(x(i,:) + k1/2, t(i) + dt/2, u_B(i,:), I_B)';
        k3 = dt * fn(x(i,:) + k2/2, t(i) + dt/2, u_B(i,:), I_B)';
        k4 = dt * fn(x(i,:) + k3, t(i) + dt, u_B(i,:), I_B)';

        x(i+1,:) = x(i,:) + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

        sigma_BN_ip1 = x(i+1,1:3);
        if norm(sigma_BN_ip1) > 1
            x(i+1,1:3) = mrp_shadow(sigma_BN_ip1);
        end
    end
    % sigma_BR = sigma_BRs;
    % omega_BR_B = omega_BRs_B;
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

