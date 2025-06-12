clear; clc; close all;

% ASEN 6060 - HW 4, Problem 1
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
mu = mass_ratio_em;

%% Part a - i)

% Get L2 Point
% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Only looking at L2 eq point planar oscillatory modes
l2_pos = [em_eq_pts(2,:), 0];

l2_in_plane_modes = in_plane_modes(mu, l2_pos);

oscillatory_eval = l2_in_plane_modes(3);
uxx_l2 = u_xx(mu, l2_pos);
uxy_l2 = u_xy(mu, l2_pos);
uyy_l2 = u_yy(mu, l2_pos);
U_star_XX = [uxx_l2, uxy_l2; uxy_l2, uyy_l2];
Omega = [0 2; -2 0];

A2D = [zeros(2), eye(2); U_star_XX, Omega];

[V, D] = eig(A2D);

oscillatory_evec = real(V(:,3));

oscillatory_pos_mag = norm([oscillatory_evec(1), oscillatory_evec(2)]);

pos_mag_req = 0.0001;

oscillatory_mag_factor = pos_mag_req / oscillatory_pos_mag;

oscillatory_ic = oscillatory_evec .* oscillatory_mag_factor;

% Time is one period
t = linspace(0, 2*pi/imag(oscillatory_eval), 1000);

xi_0 = oscillatory_ic(1);
xi_dot_0 = oscillatory_ic(3);
eta_0 = oscillatory_ic(2);
eta_dot_0 = oscillatory_ic(4);
x0 = [l2_pos(1) + xi_0; l2_pos(2) + eta_0; 0; xi_dot_0; eta_dot_0; 0];

TOL = 5e-14;

% Set options for ode113
options = odeset('RelTol', TOL, 'AbsTol', TOL);
[tout, xout] = ode113(@(t, state)CR3BP(state, mu), [0 t(end)], x0, options);

V0 = [x0; t(end)];

V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, true);

[tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

figure()
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
hold on
scatter(V_soln(1), V_soln(2), 'filled', 'blue')
plot(xout_corrected(:,1), xout_corrected(:,2), 'LineWidth',2)
hold off
axis equal
legend("L2", "Initial State", "Corrected Trajectory")
grid on
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
hold off
title("L2 Lyapunov Orbit - Corrected Trajectory")

figure()
plot(xout(:,1), xout(:,2), 'LineWidth',2)
hold on
scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
scatter(xi_0 + l2_pos(1), eta_0 + l2_pos(2), 'filled', 'blue')

plot(xout_corrected(:,1), xout_corrected(:,2), 'LineWidth',2)
axis equal
legend("CR3BP Propagation", "L2", "Initial State", "Corrected Trajectory")
grid on
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
hold off
title("L2 Lyapunov Orbit - Initial and Corrected Trajectory")

% V_family_a1 = pseudo_arc_length_continuation(V_soln, mu);
% save("V_family_a1.mat", "V_family_a1")

load("V_family_a1.mat")

p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

figure()
scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
hold on
scatter(p1_pos(1), p1_pos(2), 'filled', 'blue')
scatter(p2_pos(1), p2_pos(2), 'filled', 'black')

for i = 1:50:size(V_family_a1, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family_a1(7,i)], V_family_a1(1:6,i), options);
    plot(xout(:,1), xout(:,2), 'LineWidth',2)
end
hold off
legend("L2", "Earth", "Moon")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
grid on
axis equal
title("L2 Lyapunov Orbit Family")

figure()
for i = 1:size(V_family_a1,2)
    x = V_family_a1(1,i);
    y = V_family_a1(2,i);
    z = V_family_a1(3,i);
    x_dot = V_family_a1(4,i);
    y_dot = V_family_a1(5,i);
    z_dot = V_family_a1(6,i);
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
    C(i) = x^2 + y^2 + z^2 + 2*(1-mu)/r1 + mu/r2 - (x_dot^2 + y_dot^2 + z_dot^2);
end
plot(V_family_a1(end,:)*t_star_em/86400, C, 'o')
hold on
scatter(V_family_a1(end,1)*t_star_em/86400, C(1), 'filled', 'red')
scatter(V_family_a1(end,end)*t_star_em/86400, C(end), 'filled', 'black')
legend("Jacobi Constant", "Start", "End")
xlabel("Period [days]")
ylabel("Jacobi Constant [-]")
grid on
title("Jacobi Constant and Period along the Family")


%% Part a - ii)

x0 = [1.180462, 0, -0.0209998, 0, -0.158363, 0]';
T = 3.411921;
V = [x0; T];

% V_family_a2 = compute_and_plot_families(V, mu, l2_pos, options);

title("L2 Halo Orbit Family")





%% Part a - iii)

x0 = [1.0301513, 0,0, 0, 0.7030025,0.1552945]';
T = 4.312367;
V = [x0; T];

V_family_a3 = compute_and_plot_families(V, mu, l2_pos, options);
title("L2 Axial Orbit Family")
save("V_family_a3.mat", "V_family_a3")

%% Part b

% L2 Halo first orbit
identity_row = reshape(eye(6), [36,1]);
[tout, xout] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_a2(7,1)], [V_family_a2(1:6,1); identity_row], options);

halo_monodromy = reshape(xout(end,7:42), [6,6])';
[V_halo, D_halo] = eig(halo_monodromy);


%% Part c

% L2 Lyapunov orbit family
[tout, xout] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_a1(7,1)], [V_family_a1(1:6,1); identity_row], options);

lyapunov_monodromy = reshape(xout(end,7:42), [6,6])';
[V_lyapunov, D_lyapunov] = eig(lyapunov_monodromy);

% L2 Axial orbit family
[tout, xout] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_a3(7,1)], [V_family_a3(1:6,1); identity_row], options);

axial_monodromy = reshape(xout(end,7:42), [6,6])';
[V_axial, D_axial] = eig(axial_monodromy);


%% Functions

function V_family = compute_and_plot_families(V0, mu, l2_pos, options)
    V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, true);

    [tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

    figure()
    scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
    hold on
    scatter(V_soln(1), V_soln(2), 'filled', 'blue')
    plot3(xout_corrected(:,1), xout_corrected(:,2), xout_corrected(:,3),'LineWidth',2)
    hold off
    axis equal
    legend("L2", "Initial State", "Corrected Trajectory")
    grid on
    xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
    hold off
    title("L2 Lyapunov Orbit - Corrected Trajectory")
    
    V_family = pseudo_arc_length_continuation(V_soln, mu);
    
    p1_pos = [-mu, 0, 0];
    p2_pos = [1-mu, 0, 0];
    
    figure()
    scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
    hold on
    scatter(p1_pos(1), p1_pos(2), 'filled', 'blue')
    scatter(p2_pos(1), p2_pos(2), 'filled', 'black')
    
    for i = 1:50:size(V_family, 2)
        [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family(7,i)], V_family(1:6,i), options);
        plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2)
    end
    hold off
    legend("L2", "Earth", "Moon")
    xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
    grid on
    axis equal
end

