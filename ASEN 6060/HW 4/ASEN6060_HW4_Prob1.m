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

V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, false);

[tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

% V_family_a1 = test_psal_cont(V_soln, mu);
% save("V_family_a1.mat", "V_family_a1")
load("V_family_a1.mat")

p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];


%% Part a - ii)

x0 = [1.180462, 0, -0.0209998, 0, -0.158363, 0]';
T = 3.411921;
V0 = [x0; T];

V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, true);

[tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

% V_family_a2 = test_psal_cont(V_soln, mu);
% save("V_family_a2.mat", "V_family_a2")
load("V_family_a2.mat")
V_family_a2_south = V_family_a2;
load("V_family_a2_pos_delta_s.mat")
V_family_a2_north = V_family_a2;

figure()
p1 = scatter3(l2_pos(1), l2_pos(2), l2_pos(3), 'filled', 'red');
hold on
p2 = scatter3(p2_pos(1), p2_pos(2), p2_pos(3), 'filled', 'black');

for i = 1:250:size(V_family_a2_south, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family_a2_south(7,i)], V_family_a2_south(1:6,i), options);
    p3 = plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2, 'Color', 'red');
end

for i = 1:250:size(V_family_a2_north, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family_a2_north(7,i)], V_family_a2_north(1:6,i), options);
    p4 = plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2, 'Color', 'blue');
end
hold off
p = [p1, p2, p3, p4];
legend(p, "L2", "Moon", "Southern Halo", "Northern Halo")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)
grid on
axis equal
title("L2 Halo Orbit Family")

%% Part a - iii)

x0 = [1.0301513, 0,0, 0, 0.7030025,0.1552945]';
T = 4.312367;
V0 = [x0; T];

V_soln = gen_3d_periodic_orbit_single_shooting(V0, mu, true);

[tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)/2], V_soln(1:6), options);

V_soln = [xout_corrected(end,:)'; V_soln(end)];

delta_s_neg = -1e-3;
delta_s_pos = 1e-3;
% V_family_a3_neg_ds = pseudo_arc_length_continuation(V_soln, mu, delta_s_neg);
% V_family_a3_pos_ds = pseudo_arc_length_continuation(V_soln, mu, delta_s_pos);
% save("V_family_a3_neg_ds.mat", "V_family_a3_neg_ds");
% save("V_family_a3_pos_ds.mat", "V_family_a3_pos_ds");
% V_family_a3_neg_ds = test_psal_cont(V_soln, mu);
V_family_a3_pos_ds = test_psal_cont(V_soln, mu);
% save("V_family_a3.mat", "V_family_a3");
% load("V_family_a3.mat")
% load("V_family_a3_neg_ds.mat")
% load("V_family_a3_pos_ds.mat")
%%
figure()
scatter3(l2_pos(1), l2_pos(2), l2_pos(3), 'filled', 'red')
hold on
scatter3(p2_pos(1), p2_pos(2), p2_pos(3), 'filled', 'black')

for i = 1:50:size(V_family_a3_neg_ds, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family_a3_neg_ds(7,i)], V_family_a3_neg_ds(1:6,i), options);
    plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2)
end

for i = 1:50:size(V_family_a3_pos_ds, 2)
    [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_family_a3_pos_ds(7,i)], V_family_a3_pos_ds(1:6,i), options);
    plot3(xout(:,1), xout(:,2), xout(:,3), 'LineWidth',2)
end
hold off
legend("L2", "Moon")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
zlabel('$$\hat{z}$$','Interpreter','Latex', 'FontSize',18)
grid on
axis equal
title("L2 Axial Orbit Family")
% save("V_family_a3.mat", "V_family_a3")

%%

% L2 Halo first orbit
identity_row = reshape(eye(6), [36,1]);
[tout, xout] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_a2(7,1)], [V_family_a2(1:6,1); identity_row], options);

halo_monodromy = reshape(xout(end,7:42), [6,6])';
[V_halo, D_halo] = eig(halo_monodromy);

[s1, s2, period] = stability_index(V_family_a2_north, mu, identity_row, options);
figure()
subplot(2,1,1)
plot(period, s1, 'LineWidth', 2)
hold on
scatter(period(1), s1(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S1")
legend("S1", "Starting Point", "Critical Values")
hold off
subplot(2,1,2)
plot(period, s2, 'LineWidth', 2)
hold on
scatter(period(1), s2(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S2")
legend("S2", "Starting Point", "Critical Values")
hold off
sgtitle("S1 and S2 as a function of Period - Halo Orbit (Northern)")

%%
% [s1, s2, period] = stability_index(V_family_a2_south, mu, identity_row, options);
figure()
subplot(2,1,1)
plot(period, s1, 'LineWidth', 2)
hold on
scatter(period(1), s1(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S1")
legend("S1", "Starting Point", "Critical Values")
hold off
subplot(2,1,2)
plot(period, s2, 'LineWidth', 2)
hold on
scatter(period(1), s2(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S2")
legend("S2", "Starting Point", "Critical Values")
hold off
sgtitle("S1 and S2 as a function of Period - Halo Orbit (Southern)")

function [s1, s2, period] = stability_index(V_family, mu, identity_row, options)
    counter = 1;
    
    for i = 1:size(V_family, 2)
        [tout, xout] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family(7,i)], [V_family(1:6,i); identity_row], options);
        monodromy_i = reshape(xout(end,7:42), [6,6])';
        [evecs, evals] = eig(monodromy_i);
        
        % Pass to function
        [s1_i, s2_i] = arrange_evals(evals, evecs);
    
        if (counter==1)
            s1(counter) = s1_i;
            s2(counter) = s2_i;
        else
            diff1 = abs(s1_i - s1(counter-1));
            diff2 = abs(s1_i - s2(counter-1));
            if diff1 < diff2
                s1(counter) = s1_i;
                s2(counter) = s2_i;
            else
                s1(counter) = s2_i;
                s2(counter) = s1_i;
            end
        end
    
            
        % for j = 1:6
        %     if (abs(evals(j)-1) < 1e-4)
        %         trivial_index(counter) = j;
        %         counter = counter + 1;
        %     end
        % end
        % 
        % if trivial_index(1) == 1
        %     s1(i) = evals(3) + evals(4);
        %     s2(i) = evals(5) + evals(6);
        % elseif trivial_index(1) == 3
        %     s1(i) = evals(1) + evals(2);
        %     s2(i) = evals(5) + evals(6);
        % elseif trivial_index(1) == 5
        %     s1(i) = evals(1) + evals(2);
        %     s2(i) = evals(3) + evals(4);
        % end
    
        % s1(i) = evals(1) + evals(2);
        % s2(i) = evals(3) + evals(4);
        period(counter) = V_family(7,i);
        counter = counter + 1;
    end
end

function [s1, s2] = arrange_evals(evals_mat, evecs)
    for j = 1:6
        evals(j) = evals_mat(j,j);
    end
    % Get first element and its reciprocal and remove from array
    eval1 = evals(1);
    over_eval1 = 1/eval1;
    evals(1) = [];
    % for j = 2:6
    %     if (abs(evals(j)-over_eval1)<1e-3)
    %         pair1_index = j;
    %     end
    % end
    diff = evals - over_eval1;
    [min_diff, pair1_index] = min(abs(diff));
    pair1 = eval1 + evals(pair1_index);

    evals(pair1_index) = [];

    eval2 = evals(1);
    over_eval2 = 1/eval2;
    evals(1) = [];
    % for j = 2:4
    %     if (abs(evals(j)-over_eval2)<1e-3)
    %         pair2_index = j;
    %     end
    % end
    diff = evals - over_eval2;
    [min_diff, pair2_index] = min(abs(diff));
    pair2 = eval2 + evals(pair2_index);

    evals(pair2_index) = [];

    pair3 = evals(1) + evals(2);

    % if (abs(pair1-2)<1e-3)
    %     s1 = pair2;
    %     s2 = pair3;
    % elseif (abs(pair2-2)<1e-3)
    %     s1 = pair1;
    %     s2 = pair3;
    % else
    %     s1 = pair1;
    %     s2 = pair2;
    % end

    pairs = [pair1, pair2, pair3];
    pairs_minus_2 = pairs - 2;
    [min_pair, min_pair_ind] = min(abs(pairs_minus_2));

    pairs(min_pair_ind) = [];
    s1 = pairs(1);
    s2 = pairs(2);

end


%%

% L2 Lyapunov orbit family
[tout, xout] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_a1(7,1)], [V_family_a1(1:6,1); identity_row], options);

lyapunov_monodromy = reshape(xout(end,7:42), [6,6])';
[V_lyapunov, D_lyapunov] = eig(lyapunov_monodromy);

% [s1, s2, period] = stability_index(V_family_a1, mu, identity_row, options);
figure()
subplot(2,1,1)
plot(period, s1, 'LineWidth', 2)
hold on
scatter(period(1), s1(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S1")
legend("S1", "Starting Point", "Critical Values")
hold off
subplot(2,1,2)
plot(period, s2, 'LineWidth', 2)
hold on
scatter(period(1), s2(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S2")
legend("S2", "Starting Point", "Critical Values")
hold off
sgtitle("S1 and S2 as a function of Period - Lyapunov Orbit Family")

%%
V_family_a3 = [V_family_a3_neg_ds, flip(V_family_a3_pos_ds)];
% V_family_a3 = V_family_a3_neg_ds;
% L2 Axial orbit family
[tout, xout] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_a3(7,1)], [V_family_a3(1:6,1); identity_row], options);

[s1, s2, period] = stability_index(V_family_a3, mu, identity_row, options);
figure()
subplot(2,1,1)
plot(period, s1, 'LineWidth', 2)
hold on
scatter(period(1), s1(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S1")
legend("S1", "Starting Point", "Critical Values")
hold off
subplot(2,1,2)
plot(period, s2, 'LineWidth', 2)
hold on
scatter(period(1), s2(1), 'filled', 'red')
yline(2, '--r')
yline(-2, '--r')
xlabel("Period [-]")
ylabel("S2")
legend("S2", "Starting Point", "Critical Values")
hold off
sgtitle("S1 and S2 as a function of Period - Axial Orbit Family")