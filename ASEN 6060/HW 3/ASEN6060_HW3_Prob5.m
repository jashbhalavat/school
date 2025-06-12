clear; clc; close all;

% ASEN 6060 - HW 3 Problem 5
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

%% Part a - Correction (from Problem 4)

% Given
x0 = [0.82340, 0, -0.026755,0,0.13742,0]';
T = 2.7477;
V0 = [x0; T];

mu = mass_ratio_em;

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Only looking at L1 eq point planar oscillatory modes
l1_pos = [em_eq_pts(1,:), 0];

% Set options for ode113
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Call ode113 function
[tout, xout] = ode113(@(t, state)CR3BP(state, mu), [0 T], x0, options);

statef_V0 = xout(end,:);

V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, false);

%% Part a - Continuation (Pseudo-arc length)

V_family = pseudo_arc_length_continuation(V_soln, mu);

%% Plot all family members

p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure()
for i = 1:size(V_family,2)
    x = V_family(1,i);
    y = V_family(2,i);
    z = V_family(3,i);
    x_dot = V_family(4,i);
    y_dot = V_family(5,i);
    z_dot = V_family(6,i);
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
    C(i) = x^2 + y^2 + z^2 + 2*(1-mu)/r1 + mu/r2 - (x_dot^2 + y_dot^2 + z_dot^2);
end
plot(V_family(end,:)*t_star_em/86400, C, 'o')
hold on
scatter(V_family(end,1)*t_star_em/86400, C(1), 'filled', 'red')
scatter(V_family(end,end)*t_star_em/86400, C(end), 'filled', 'black')
legend("Jacobi Constant", "Start", "End")
xlabel("Period [days]")
ylabel("Jacobi Constant [-]")
grid on
title("Jacobi Constant and Period along the Family")

figure()
subplot(1,2,1)
scatter3(l1_pos(1), l1_pos(2), l1_pos(3), 'filled', 'red')
hold on
scatter3(p2_pos(1), p2_pos(2), p2_pos(3), 'filled', 'black')

for i = 1:100:size(V_family, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family(7,i)], V_family(1:6,i), options);
    plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2)
end
hold off
legend("L1", "Moon")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
title("L1 Halo Orbit Family")
axis equal

subplot(1,2,2)
scatter3(l1_pos(1), l1_pos(2), l1_pos(3), 'filled', 'red')
hold on
scatter3(p2_pos(1), p2_pos(2), p2_pos(3), 'filled', 'black')

for i = 1:100:size(V_family, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family(7,i)], V_family(1:6,i), options);
    plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2)
end
hold off
legend("L1", "Moon")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)
title("L1 Halo Orbit Family")
axis equal
