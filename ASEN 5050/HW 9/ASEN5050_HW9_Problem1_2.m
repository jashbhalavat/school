clear; clc; close all;

% ASEN 5050 - HW 9 - Problem 1, 2
% Fall 2024
% Jash Bhalavat

%% Given
a_sc = 10000; % km
e = 0;
mu_earth = 3.986004415e5; % km^3/s^2
r_earth = 6378.1363; % km

P = 2*pi*sqrt(a_sc^3/mu_earth);

%% Problem 1

% Cubesat 1 - prox ops
cs1_x0 = 5; % m
cs1_xdot0 = 0; % m/s

n = sqrt(mu_earth/a_sc^3); % rad/s

cs1_y0 = 5; % m
cs1_ydot0 = -2*n*cs1_x0; % m/s

t = 0:10:(2*pi/n);

tspan = [0,pi/(2*n), pi/n, 3*pi/(2*n), 2*pi/n];

for i = 1:length(t)
    x_t(i) = cs1_x0 * cos(n*t(i));
    y_t(i) = cs1_y0 - 2*cs1_x0*sin(n*t(i));
end

% plot(y_t, x_t)

%% Problem 2

% @t = 0
cs2_x0 = 2; % m
cs2_y0 = 2; % m 
cs2_z0 = 0; % m 
cs2_xdot0 = -0.03; % m/s
cs2_ydot0 = 0.01; % m/s 
cs2_zdot0 = 0.05; % m/s 

% @t1 = P/2
t1 = P/2;
cs2_x1 = -2; % m
cs2_y1 = 2; % m 
cs2_z1 = 0; % m 
cs2_xdot1 = 0; % m/s
cs2_ydot1 = 0; % m/s 
cs2_zdot1 = 0; % m/s 

% r0 = [cs2_x0; cs2_y0; cs2_z0];
r0 = [cs2_x0; cs2_y0];
% rf = [cs2_x1; cs2_y1; cs2_z1];
rf = [cs2_x1; cs2_y1];
cs2_v0_plus = phi_rv(n, t1) \ (rf - phi_rr(n,t1)*r0);
cs2_v0_plus(3) = 0;
delta_v_1 = cs2_v0_plus - [cs2_xdot0, cs2_ydot0, cs2_zdot0]';
norm_dv1 = norm(delta_v_1);

test = [r0;0;cs2_v0_plus;0];
post_dv1 = cw(test, n, t1)

% Part b
cs2_v1_minus = phi_vr(n,t1)*r0 + phi_vv(n,t1)*cs2_v0_plus(1:2);
cs2_v1_minus(3) = 0;
delta_v_2 = cs2_v1_minus - [cs2_xdot1, cs2_ydot1, cs2_zdot1]';
norm_dv2 = norm(delta_v_2);

% Part c
tvec = 0:10:t1;
cs2_transfer_x0 = [cs2_x0, cs2_y0, cs2_z0, cs2_v0_plus(1), cs2_v0_plus(2), cs2_v0_plus(3)];
for i = 1:length(tvec)
    cs2_transfer_state(:,i) = cw(cs2_transfer_x0, n, tvec(i));
end

plot(tvec, cs2_transfer_state)


function out = phi_rr(n, t)
    out = zeros(2,2);
    out(1,1) = 4 - 3*cos(n*t);
    out(2,1) = 6 * (sin(n*t) - n*t);
    out(2,2) = 1;
    % out(3,3) = cos(n*t);
end

function out = phi_rv(n, t)
    out = zeros(2,2);
    out(1,1) = 1/n * sin(n*t);
    out(2,1) = 2/n * (cos(n*t) - 1);
    out(1,2) = 2/n * (1 - cos(n*t));
    out(2,2) = 4/n * sin(n*t) - 3*t;
    % out(3,3) = 1/n * sin(n*t);
end

function out = phi_vr(n, t)
    out = zeros(2,2);
    out(1,1) = 3*n*sin(n*t);
    out(2,1) = 6*n*(cos(n*t)-1);
end

function out = phi_vv(n, t)
    out = zeros(2,2);
    out(1,1) = cos(n*t);
    out(2,1) = -2*sin(n*t);
    out(1,2) = 2*sin(n*t);
    out(2,2) = 4*cos(n*t) - 3;
end

function out = cw(x0, n, t)
    out(1) = 4*x0(1) + 2/n*x0(5) + x0(4)/n*sin(n*t) - (3*x0(1) + 2/n*x0(5)) * cos(n*t);
    out(2) = x0(2) - 2/n*x0(4) - 3*(2*n*x0(1) + x0(5))*t + 2 * (3*x0(1) + 2/n*x0(5)) * sin(n*t) + 2/n*x0(4)*cos(n*t);
    out(3) = 1/n * x0(6) * sin(n*t) + x0(3)*cos(n*t);
end



