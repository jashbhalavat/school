clear; clc; close all;

% ASEN 5044 - Final Project - Part 1
% Jash Bhalavat
% Fall 2024

r0 = 6678; % km
mu = 398600; % km^3/s^2
k = sqrt(mu/r0^3); % mean motion

n = 4;
m = 2;
p = 3;

%% Question 1
% A_nom = [0 1 0 0; (2*mu/r0^3) 0 0 0; 0 0 0 1; 0 0 -mu/r0^3 0];
% B_nom = [0 0; 1 0; 0 0; 0 1];

%% Question 2
dt = 10; % s

% Fh = expm([A_nom, B_nom; zeros(2, 6)]*dt);
% 
% F = Fh(1:4, 1:4);
% G = Fh(1:4, 5:6);

%% Question 3

% Time steps
t = 0:dt:14000;
dx = zeros(4, length(t));
u = zeros(2, length(t));

xnom_0 = [r0; 0; 0; r0*sqrt(mu/r0^3)];
x_nom_k = [r0 * cos(k.*t);
    -r0 * k * sin(k.*t);
    r0 * sin(k.*t);
    r0 * k * cos(k.*t)];

% From progress report #1
dx0 = [0; 0.075; 0; -0.021];
dx = [dx0];

I_n = eye(4);

for i = 2:length(t)
    % x1_k_minus_1 = x_nom_k(1,i-1);
    % x3_k_minus_1 = x_nom_k(3,i-1);
    % x1_k_minus_1 = xnom_0(1);
    % x3_k_minus_1 = xnom_0(3);
    % t_i = t(i);
    % t_i_minus_1 = t(i-1);
    % x_nom_k(:,i) = [-mu*x1_k_minus_1*t_i^2/(2*r0^3) + r0;
    %     -mu*x1_k_minus_1*t_i/(r0^3);
    %     -mu*x3_k_minus_1*t_i^2/(2*r0^3) + r0*sqrt(mu/(r0^3))*t_i;
    %     -mu*x3_k_minus_1*t_i/(r0^3) + r0*sqrt(mu/r0^3)];
    
    F_k = I_n;
    dx(:,i) = F_k * dx(:,i-1);
end


x = x_nom_k + dx;



% % From progress report
% dx0 = [0; 0.075; 0; -0.021];
% 
% % Initialize perturbation state
% dx(:, 1) = dx0;
% % Compute perturbation states over time
% for i = 2:length(t)
%     dx(:, i) = F*dx(:, i - 1) + G*u(:, i);
% end

figure();
subplot(4, 1, 1);
plot(t, x(1, :));
title('Linearized approx perturbations vs. Time');
ylabel('\deltaX (km)');
grid on;
subplot(4, 1, 2);
plot(t, x(2, :));
ylabel('\deltaX-dot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, x(3, :));
ylabel('\deltaY (km)');
grid on;
subplot(4, 1, 4);
plot(t, x(4, :));
ylabel('Y-dot (km/s)');
xlabel('\deltaTime t (s)');
grid on;

%% ODE45 NL Dynamics
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


