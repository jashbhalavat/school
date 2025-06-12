clear; clc; close all;

% ASEN 5010 - HW 6, CC1 (sat ctrl), Problem 4
% Spring 2025
% Jash Bhalavat

I_B = diag([100, 75, 80]); % kg-m2
sigma_BN_t0 = [0.1, 0.2, -0.1]';
omega_BN_B_t0 = [30, 10, -20]' * pi/180;
K = 5;
P = 10 * eye(3);
t_span = 0:0.01:120.02;
f = 0.05; % rad/s

state_0 = [sigma_BN_t0; omega_BN_B_t0];

L = [0.0, 0.0, 0.0]';

[state, u, sigma_BR] = rk4(@(state, t, u_B, I_B)diff_eq(state, t, u_B, I_B, L), t_span, state_0, I_B, K, P, f, L);

sigma_BR_at_60 = sigma_BR(:, 6001);
sigma_BR_norm_at_60 = sqrt(sigma_BR_at_60(1)^2 + sigma_BR_at_60(2)^2 + sigma_BR_at_60(3)^2);

function state_dot = diff_eq(state, t, u, I, L)
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

    omega_BN_B_dot(1,1) = -(I(3,3) - I(2,2))/I(1,1) * omega_BN_B(2) * omega_BN_B(3) + u(1)/I(1,1) + L(1)/I(1,1);
    omega_BN_B_dot(2,1) = -(I(1,1) - I(3,3))/I(2,2) * omega_BN_B(1) * omega_BN_B(3) + u(2)/I(2,2) + L(2)/I(2,2);
    omega_BN_B_dot(3,1) = -(I(2,2) - I(1,1))/I(3,3) * omega_BN_B(2) * omega_BN_B(1) + u(3)/I(3,3) + L(3)/I(3,3);
    
    state_dot = [sigma_BN_dot; omega_BN_B_dot];
end

function [x, u_B, sigma_BR] = rk4(fn, t, x0, I_B, K, P, f, L)

    dt = t(2) - t(1);
    
    % x is nx6
    x = x0';

    for i = 1:length(t)-2
        sigma_RN = [0.2*sin(f*t(i)); 0.3*cos(f*t(i)); -0.3*sin(f*t(i))];
        RN = mrp_to_dcm(sigma_RN);

        sigma_RN_ip1 = [0.2*sin(f*t(i+1)); 0.3*cos(f*t(i+1)); -0.3*sin(f*t(i+1))];
        RN_ip1 = mrp_to_dcm(sigma_RN_ip1);

        sigma_RN_ip2 = [0.2*sin(f*t(i+2)); 0.3*cos(f*t(i+2)); -0.3*sin(f*t(i+2))];
        RN_ip2 = mrp_to_dcm(sigma_RN_ip2);

        RN_dot = (RN_ip1 - RN) / dt;
        omega_RN_R_tilde = -RN_dot * RN';
        omega_RN_R = [-omega_RN_R_tilde(2,3); omega_RN_R_tilde(1,3); -omega_RN_R_tilde(1,2)];

        RN_dot_ip1 = (RN_ip2 - RN_ip1) / dt;
        omega_RN_R_tilde_ip1 = -RN_dot_ip1 * RN_ip1';
        omega_RN_R_ip1 = [-omega_RN_R_tilde_ip1(2,3); omega_RN_R_tilde_ip1(1,3); -omega_RN_R_tilde_ip1(1,2)];

        sigma_BN = x(i,1:3)';
        BN = mrp_to_dcm(sigma_BN);
        BR = BN * RN';
        sigma_BR(:,i) = dcm_to_mrp(BR);
        
        omega_BN_B = x(i,4:6)';
        omega_RN_B = BR * omega_RN_R;
        omega_BR_B = omega_BN_B - omega_RN_B;

        omega_dot_RN_R = (omega_RN_R_ip1 - omega_RN_R) / dt;
        omega_dot_RN_B = BR * omega_dot_RN_R;
        
        u_B(i,:) = (-K*sigma_BR(:,i) -P*omega_BR_B + I_B*(omega_dot_RN_B - cross(omega_BN_B, omega_RN_B)) + skew_symmetric(omega_BN_B)*I_B*omega_BN_B)';

        for j = 1:3
            if abs(u_B(i,j)) >= 1
                u_B(i,j) = sign(u_B(i,j));
            end
        end

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

