clear; clc; close all;

% ASEN 6060 - Week 10 Discussion
% Spring 2025
% Jash Bhalavat

%% Constants

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

global count poincare_stored

%%

mu = mass_ratio_em;

x0_1 = [0.3693,0,0,0,1.4772,0];
x0_2 = [0.3670,0,0,0,1.4865,0];
x0_3 = [0.3640,0,0,0,1.4994,0];

t1 = 5000;
t2 = 5000;
t3 = 5000;

n_crossings = 2000;

TOL = 5e-14;

% Set options for ode113
options = odeset('RelTol', TOL, 'AbsTol', TOL, 'Events', @(t,y)eventFn(t,y,mu,n_crossings));

count = 0;
poincare_stored = [];
[tout_1, xout_1] = ode113(@(t,state)CR3BP(state, mu), [0, t1], x0_1, options);
figure()
scatter(poincare_stored(:,1), poincare_stored(:,2), 3, 'filled', 'black');
ylabel("x")
xlabel("$\dot{x}$", 'Interpreter','latex')
title("Poincar\'e Map", 'Interpreter','latex')

count = 0;
poincare_stored = [];
[tout_2, xout_2] = ode113(@(t,state)CR3BP(state, mu), [0, t2], x0_2, options);
[tout_1, xout_1] = ode113(@(t,state)CR3BP(state, mu), [0, t1], x0_1, options);
figure()
scatter(poincare_stored(:,1), poincare_stored(:,2), 3, 'filled', 'black');
ylabel("x")
xlabel("$\dot{x}$", 'Interpreter','latex')
title("Poincar\'e Map", 'Interpreter','latex')

count = 0;
poincare_stored = [];
[tout_3, xout_3] = ode113(@(t,state)CR3BP(state, mu), [0, t3], x0_3, options);
[tout_1, xout_1] = ode113(@(t,state)CR3BP(state, mu), [0, t1], x0_1, options);
figure()
scatter(poincare_stored(:,1), poincare_stored(:,2), 3, 'filled', 'black');
ylabel("x")
xlabel("$\dot{x}$", 'Interpreter','latex')
title("Poincar\'e Map", 'Interpreter','latex')

% figure()
% plot(xout_1(:,1), xout_1(:,2), 'LineWidth', 2)
% hold on
% scatter(-mu, 0, 'filled', 'blue')
% scatter(1-mu, 0, 'filled', 'black')
% hold off
% grid on
% xlabel("x")
% ylabel("y")
% title("x_{0,1}")
% 
% figure()
% plot(xout_2(:,1), xout_2(:,2), 'LineWidth', 2)
% hold on
% scatter(-mu, 0, 'filled', 'blue')
% scatter(1-mu, 0, 'filled', 'black')
% hold off
% grid on
% xlabel("x")
% ylabel("y")
% title("x_{0,2}")
% 
% figure()
% plot(xout_3(:,1), xout_3(:,2), 'LineWidth', 2)
% hold on
% scatter(-mu, 0, 'filled', 'blue')
% scatter(1-mu, 0, 'filled', 'black')
% hold off
% grid on
% xlabel("x")
% ylabel("y")
% title("x_{0,3}")

function state_dot = CR3BP(state, mu)
    % Circular Restricted 3 Body Problem non-dimensional EOMs
    x = state(1);
    y = state(2);
    z = state(3);
    xdot = state(4);
    ydot = state(5);
    zdot = state(6);

    r1 = sqrt((x + mu)^2 + (y)^2 + (z)^2);
    r2 = sqrt((x - 1 + mu)^2 + (y)^2 + (z)^2);

    state_dot(1, 1) = xdot;
    state_dot(2, 1) = ydot;
    state_dot(3, 1) = zdot;

    state_dot(4, 1) = 2*ydot + x - (1 - mu)*(x + mu)/(r1^3) - mu * (x - 1 + mu)/(r2^3);
    state_dot(5, 1) = -2*xdot + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3);
    state_dot(6, 1) = - (1 - mu)*z/(r1^3) - mu*z/(r2^3);
end

function [value,isterminal,direction] = eventFn(t,y,mu,n_crossings)
    global count;
    global poincare_stored;
    if count < n_crossings
        value = y(2);
        isterminal = 0;
        direction = 1;
        if (abs(value) < 5e-14 && y(5) > 0)
            count = count + 1;
            poincare_stored = [poincare_stored; y(1), y(4)];
        end
    elseif count == n_crossings
        value = y(2); % Want x to be 1-mu
        isterminal = 1; % Halt integration when value is 0
        direction = 1; % When zero is approached from +ve i.e. x_dot > 0
    end
end
