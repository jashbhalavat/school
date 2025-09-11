clear; clc; close all;

% ASEN 6014 - Fall 2025
% HW 2, Quiz 5, Problem 2
% CWH Equations
% Jash Bhalavat

%% Given

mu = 398600; % km3/s2
a = 7500; % km

yoff = 200; % km

r_c_H = [a, 0, 0]'; % km
v_c_H = [0, sqrt(mu/a), 0]'; % km/s

rho_H_0 = [0, yoff, 0]'; % km
rho_dot_H_0 = [0, 0, 0]'; % km

% Assume H = N for chief
r_c_N = r_c_H;
v_c_N = v_c_H;

state_0 = [r_c_N; v_c_N; rho_H_0; rho_dot_H_0];

time = 2000; % sec

[tout, state_out] = ode45(@(t,state)nonlinear_rel_eom(t,state), [0, time], state_0);

function [r_d_N, v_d_N] = hill_rel_2_inertial_deputy(r_c_N, v_c_N, rho_O, rho_O_prime)
   ohat_r = r_c_N/norm(r_c_N);
   h_N = cross(r_c_N, v_c_N);
   ohat_h = h_N/norm(h_N);
   ohat_theta = cross(ohat_h, ohat_r);

   dcm_ON = [ohat_r'; ohat_theta'; ohat_h'];
   dcm_NO = dcm_ON';

   fdot = norm(h_N)/(norm(r_c_N))^2;
   omega_ON_O = [0; 0; fdot];

   r_d_N = r_c_N + dcm_NO*rho_O;
   v_d_N = v_c_N + dcm_NO * (rho_O_prime + cross(omega_ON_O, rho_O));
end

sep_dist_tf = norm(state_out(end, 7:9)) - yoff;

function state_dot = nonlinear_rel_eom(t, state)
    % State - [r_c_N, v_c_N, rho_H, rho_dot_H]

    % chief EOMs
    mu = 398600; % km3/s2
    r_c_N = state(1:3);
    v_c_N = state(4:6);
    rho_O = state(7:9);
    rho_O_prime = state(10:12);
    rc = norm(r_c_N);
    ohat_r = r_c_N/norm(r_c_N);
    rc_dot = dot(v_c_N, ohat_r); % ASSUMPTION!!!
    r_c_N_dotdot = -mu/rc^3 * r_c_N;
    state_dot(1:3, 1) = v_c_N;
    state_dot(4:6, 1) = r_c_N_dotdot;

    [r_d_N, v_d_N] = hill_rel_2_inertial_deputy(r_c_N, v_c_N, rho_O, rho_O_prime);
    rd = norm(r_d_N);

    % rho EOMs
    x = state(7);
    y = state(8);
    z = state(9);
    xdot = state(10);
    ydot = state(11);
    zdot = state(12);

    h = norm(cross(r_c_N, v_c_N));
    fdot = h/rc^2;
    x_dotdot = -mu/rd^3 * (rc + x) + 2*fdot*(ydot - y*rc_dot/rc) + x*fdot^2 + mu/rc^2;
    y_dotdot = -mu/rd^3*y - 2*fdot*(xdot - x*rc_dot/rc) + y*fdot^2;
    z_dotdot = -mu/rd^3*z;

    state_dot(7:9,1) = [xdot; ydot; zdot];
    state_dot(10:12,1) = [x_dotdot; y_dotdot; z_dotdot];

end
