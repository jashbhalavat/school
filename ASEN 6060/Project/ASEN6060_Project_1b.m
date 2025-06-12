clear; clc; close all;

% ASEN 6060 - Final Project, 1b
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

p1_pos = [-mu, 0, 0];
p2_pos = [1-mu, 0, 0];

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Only looking at L2 eq point planar oscillatory modes
l1_pos = [em_eq_pts(1,:), 0];
l2_pos = [em_eq_pts(2,:), 0];

TOL = 5e-14;

% Set options for ode113
options = odeset('RelTol', TOL, 'AbsTol', TOL);

global count poincare_stored counted_once

%% 1b Orbits

% L1 orbit is a Lyapunov orbit
load("V_family_L1_Lyapunov.mat")
V_family_L1_Lyapunov = V_family;
l1_orbit_idx = 50;

% L2 orbit is a Lyapunov orbit
load("V_family_L2_Lyapunov.mat")
V_family_L2_Lyapunov = V_family_a1;
l2_orbit_idx = 10;

figure()
scatter(l1_pos(1), l1_pos(2), 'filled', 'blue')
hold on
scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
scatter(p2_pos(1), p2_pos(2), 'filled', 'black')

[tout_L1, xout_L1] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_L1_Lyapunov(7,l1_orbit_idx)], [V_family_L1_Lyapunov(1:6,l1_orbit_idx); reshape(eye(6), [36,1])], options);
plot(xout_L1(:,1), xout_L1(:,2), 'LineWidth',2)

[tout_L2, xout_L2] = ode113(@(t,state)CR3BP_full(state, mu), [0, V_family_L2_Lyapunov(7,l2_orbit_idx)], [V_family_L2_Lyapunov(1:6,l2_orbit_idx); reshape(eye(6), [36,1])], options);
plot(xout_L2(:,1), xout_L2(:,2), 'LineWidth',2)

cj_L1 = jacobi_constant(V_family_L1_Lyapunov(:,l1_orbit_idx), mu);
cj_L2 = jacobi_constant(V_family_L2_Lyapunov(:,l2_orbit_idx), mu);
delta_cj = cj_L2 - cj_L1;

hold off
legend("L1", "L2", "Moon", "L1 Lyapunov Orbit", "L2 Lyapunov Orbit")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
grid on
axis equal
title("1b Initial & Final Periodic Orbits")

function out = u_star_times_2(x, y, z, mu)
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    out = (x^2 + y^2) + 2*(1 - mu)/r1 + 2*mu/r2;
end

function C = jacobi_constant(x, mu)
    C = u_star_times_2(x(1), x(2), x(3), mu) - x(4)^2 - x(5)^2 - x(6)^2;
end


%% Manifolds 
n_crossings = 2;
L1_manifold_time = 10;
L2_manifold_time = 6;

L1_manifolds = part_b(tout_L1, xout_L1, mu, l1_pos, L1_manifold_time, n_crossings);

poincare_stored = [];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) b_eventFn(t, y, mu, n_crossings));

plot(xout_L1(:,1), xout_L1(:,2), 'black', 'LineWidth', 3)
hold on
plot(xout_L2(:,1), xout_L2(:,2), 'black', 'LineWidth', 3)
scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
scatter(p2_pos(1), p2_pos(2), 'filled', ' black')

for k = 1:length(L1_manifolds)
    count = 0;
    counted_once = false;
    [moon_stable_t, moon_stable_x] = ode113(@(t, state)CR3BP(state, mu), [0, -L1_manifold_time], L1_manifolds(:,k), options);
    stable_end_times(k) = moon_stable_t(end);
    stable_end_states(k, :) = moon_stable_x(end, :);
    plot(moon_stable_x(:,1), moon_stable_x(:,2), 'blue')
end
grid on
axis equal
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)    

poincare_stable = poincare_stored;
%%

L2_manifolds = part_c(tout_L2, xout_L2, mu, l2_pos, L2_manifold_time, n_crossings);

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) c_eventFn(t, y, mu, n_crossings));
poincare_stored = [];

for k = 1:length(L2_manifolds)
    count = 0;
    counted_once = false;
    [moon_unstable_t, moon_unstable_x] = ode113(@(t, state)CR3BP(state, mu), [0, L2_manifold_time], L2_manifolds(:,k), options);
    unstable_end_times(k) = moon_unstable_t(end);
    unstable_end_states(k, :) = moon_unstable_x(end, :);
    plot(moon_unstable_x(:,1), moon_unstable_x(:,2), 'red')
end
poincare_unstable = poincare_stored;

hold off
title("Moon-Bound Stable/Unstable Manifolds associated with L1, L2 Lyapunov Orbits")


%% Poincare Map

figure()
scatter(poincare_unstable(:,2), poincare_unstable(:,1), 10, 'filled', 'red');
hold on
scatter(poincare_stable(:,2), poincare_stable(:,1), 10, 'filled', 'blue');
xlabel("$\dot{y}$", 'Interpreter','latex')
ylabel("y")
title("Poincar\'e Map", 'Interpreter','latex')
grid on
hold off
legend("Unstable", "Stable")

%% Multiple Shooting Initial Guess

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

unstable_ig = 127;
stable_ig = 89;

% Initial guess
x_1_0 = xout_L2(1, 1:6)';
delta_T1 = tout_L2(unstable_ig);
% [tout, x_1] = ode113(@(t, state)CR3BP(state, mu), [0, delta_T1], x_1_0);
x_1_f = xout_L2(unstable_ig, 1:6)';

x_2_0 = L2_manifolds(:, unstable_ig);
delta_T2 = stable_end_times(unstable_ig);
x_2_f = stable_end_states(unstable_ig, :)';

x_3_f = L1_manifolds(:, stable_ig);
delta_T3 = -stable_end_times(stable_ig);
x_3_0 = stable_end_states(stable_ig, :)';

x_4_0 = xout_L2(stable_ig, 1:6)';
delta_T4 = V_family_L2_Lyapunov(7, l2_orbit_idx) - tout_L2(stable_ig);
[tout, x_4] = ode113(@(t, state)CR3BP(state, mu), [0, delta_T4], x_4_0);
x_4_f = x_4(end, 1:6)';

r_des_4 = xout_L2(1, 1:3)';

V0 = [x_1_0; delta_T1; x_2_0; delta_T2; x_3_0; delta_T3; x_4_0; delta_T4];

e1 = x_1_f - x_1_0;
e2 = x_2_f - x_2_0;
e3 = x_3_f - x_3_0;
F = [e1(1:3); e2(1:3); e3(1:3); x_4_f(1:3) - r_des_4];

V_soln = multiple_shooting(V0, mu, true, r_des_4);

% x_1_0 = xout_L1(unstable_ig, 1:6)';
% ig_um = L1_manifolds(:, unstable_ig);
% ig_sm = L2_manifolds(:, stable_ig);
% ig_L2_PO = xout_L2(stable_ig, 1:6)';

%% Uncorrected Trajectory

[~, x_1] = ode113(@(t, state)CR3BP(state, mu), [0, V0(7)], V0(1:6), options);
[~, x_2] = ode113(@(t, state)CR3BP(state, mu), [0, V0(14)], V0(8:13), options);
[~, x_3] = ode113(@(t, state)CR3BP(state, mu), [0, V0(21)], V0(15:20), options);
[~, x_4] = ode113(@(t, state)CR3BP(state, mu), [0, V0(28)], V0(22:27), options);

figure()
scatter(l1_pos(1), l1_pos(2), 'filled', 'blue')
hold on
scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
scatter(p2_pos(1), p2_pos(2), 'filled', 'black')

% plot(x_1(:,1), x_1(:,2), 'LineWidth', 2)
plot(x_2(:,1), x_2(:,2), 'LineWidth', 2)
% plot(x_3(:,1), x_3(:,2), 'LineWidth', 2)
% plot(x_4(:,1), x_4(:,2), 'LineWidth', 2)

hold off
legend("L1", "L2", "Moon")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
grid on
axis equal
title("Uncorrected Trajectories")

%% Corrected Trajectory

[~, x_1] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(7)], V_soln(1:6), options);
[~, x_2] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(14)], V_soln(8:13), options);
[~, x_3] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(21)], V_soln(15:20), options);
[~, x_4] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(28)], V_soln(22:27), options);

x1f = x_1(end,:);

figure()
scatter(l1_pos(1), l1_pos(2), 'filled', 'blue')
hold on
scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
scatter(p2_pos(1), p2_pos(2), 'filled', 'black')

plot(x_1(:,1), x_1(:,2), 'LineWidth', 2)
plot(x_2(:,1), x_2(:,2), 'LineWidth', 2)
plot(x_3(:,1), x_3(:,2), 'LineWidth', 2)
plot(x_4(:,1), x_4(:,2), 'LineWidth', 2)

hold off
legend("L1", "L2", "Moon")
xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
grid on
axis equal
title("Corrected Trajectories")


%% Functions 


function moon_bound_stable = part_b(tout, xout, mu, l1_pos, manifold_time, n_crossings)
    % Set options for ode113()
    % Part b
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) eventFn(t, y, mu));
    
    a = 384400; % [kg] EM average SMA
    d = 50 / a; % [-] Unitless, normalized by a
    
    period = tout(end);

    p1_pos = [-mu, 0, 0];
    p2_pos = [1-mu, 0, 0];
    
    % figure()
    % plot(xout(:,1), xout(:,2), 'black', 'LineWidth', 3)
    % hold on
    % scatter(l1_pos(1), l1_pos(2), 'filled', 'red')
    % scatter(p1_pos(1), p1_pos(2), 'filled', 'blue')
    % scatter(p2_pos(1), p2_pos(2), 'filled', ' black')

    % Compute STM - phi(t1+T, t1)
    phi_t1T_t1 = reshape(xout(end,7:42), [6,6])';

    moon_stable_cnt = 0;

    % Begin for loop
    for i = 1:length(tout)
        
        % Compute STM - phi(tj+T, tj)
        phi_tj_t1 = reshape(xout(i, 7:42), [6,6])';
        phi_tjT_tj = phi_tj_t1 * phi_t1T_t1 * inv(phi_tj_t1);

        % Get evals, evecs
        [V, D] = eig(phi_tjT_tj);

        % Get evals as an array
        for j = 1:6
            evals(j) = D(j,j);
        end

        % Subtract evals by 1 and get 2 minimum indices. These are trivial
        % indices
        evals_minus_1 = evals - 1;
        [min_evals, trivial_index] = mink(abs(evals_minus_1), 2);

        % If eval is real and not trivial, assign stable and unstable
        % indices
        for j = 1:2
            if (isreal(evals(j)) && isnotin(trivial_index, j))
                if evals(j) < 1
                    stable_index = j;
                elseif evals(j) > 1
                    unstable_index = j;
                end
            end
        end

        % Get unstable evec and normalize eigenvector by 1st 3 terms
        stable_eval = D(stable_index, stable_index);
        stable_evec = V(:, stable_index);
        stable_pos_norm = norm(stable_evec(1:3));
        stable_evec = stable_evec/stable_pos_norm;

        % ONLY FOR L1
        % If x-velocity is positive, moon-bound
        % If x-velocity if negative, earth-bound
        x_manifold_s_p = xout(i,1:6)' + d * stable_evec;
        x_manifold_s_n = xout(i,1:6)' - d * stable_evec;
        if (x_manifold_s_p(4) > 0)
            moon_stable = x_manifold_s_p;
            earth_stable = x_manifold_s_n;
        else
            moon_stable = x_manifold_s_n;
            earth_stable = x_manifold_s_p;
        end
        
        % Propagate using the event functions
        [moon_stable_t, moon_stable_x] = ode113(@(t, state)CR3BP(state, mu), [0, -manifold_time], moon_stable, options);
        [earth_stable_t, earth_stable_x] = ode113(@(t, state)CR3BP(state, mu), [0, -manifold_time], earth_stable, options);

        % plot(moon_unstable_x(:,1), moon_unstable_x(:,2), 'red')
        % plot(earth_unstable_x(:,1), earth_unstable_x(:,2), 'red')

        if abs(moon_stable_x(end,1) - (1-mu)) < 1e-6
            moon_stable_cnt = moon_stable_cnt + 1;
            moon_bound_stable(:,moon_stable_cnt) = moon_stable;
        elseif abs(earth_stable_x(end,1) - (1-mu)) < 1e-6
            moon_stable_cnt = moon_stable_cnt + 1;
            moon_bound_stable(:,moon_stable_cnt) = earth_stable;
        end

    end

    % global count;
    % global poincare_stored;
    % poincare_stored = [];
    % for k = 1:moon_unstable_cnt
    %     count = 0;
    %     options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) b_eventFn(t, y, mu, n_crossings));
    %     [moon_unstable_t, moon_unstable_x] = ode113(@(t, state)CR3BP(state, mu), [0, manifold_time], moon_bound_unstable(:,k), options);
    %     plot(moon_unstable_x(:,1), moon_unstable_x(:,2), 'red')
    % end
    % hold off
    % legend("Lyapunov Orbit", "L1", "Earth", "Moon")
    % grid on
    % axis equal
    % xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    % ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)    
end

function moon_bound_stable = part_c(tout, xout, mu, l2_pos, manifold_time, n_crossings)
    % Set options for ode113()
    % Part c
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) eventFn(t, y, mu));
    
    a = 384400; % [kg] EM average SMA
    d = 50 / a; % [-] Unitless, normalized by a
    
    period = tout(end);

    p1_pos = [-mu, 0, 0];
    p2_pos = [1-mu, 0, 0];
    
    % figure()
    % plot(xout(:,1), xout(:,2), 'black', 'LineWidth', 3)
    % hold on
    % scatter(l2_pos(1), l2_pos(2), 'filled', 'red')
    % scatter(p1_pos(1), p1_pos(2), 'filled', 'blue')
    % scatter(p2_pos(1), p2_pos(2), 'filled', ' black')

    % Compute STM - phi(t1+T, t1)
    phi_t1T_t1 = reshape(xout(end,7:42), [6,6])';

    moon_unstable_cnt = 0;

    % Begin for loop
    for i = 1:length(tout)
        % print_out = sprintf('Discrete Time Step  - %d', i);
        % disp(print_out)
        
        % Compute STM - phi(tj+T, tj)
        phi_tj_t1 = reshape(xout(i, 7:42), [6,6])';
        phi_tjT_tj = phi_tj_t1 * phi_t1T_t1 * inv(phi_tj_t1);

        % Get evals, evecs
        [V, D] = eig(phi_tjT_tj);

        % Get evals as an array
        for j = 1:6
            evals(j) = D(j,j);
        end

        % Subtract evals by 1 and get 2 minimum indices. These are trivial
        % indices
        evals_minus_1 = evals - 1;
        [min_evals, trivial_index] = mink(abs(evals_minus_1), 2);

        % If eval is real and not trivial, assign stable and unstable
        % indices
        % for j = 1:6
        %     if (isreal(evals(j)) && isnotin(trivial_index, j))
        %         if evals(j) < 1
        %             stable_index = j;
        %         elseif evals(j) > 1
        %             unstable_index = j;
        %         end
        %     end
        % end
        [~, unstable_index] = max(real(evals_minus_1));

        % Get stable/unstable evec and normalize eigenvector by 1st 3 terms
        % stable_eval = D(stable_index, stable_index);
        unstable_evec = V(:, unstable_index);
        unstable_pos_norm = norm(unstable_evec(1:3));
        unstable_evec = unstable_evec/unstable_pos_norm;
        % stable_evec(4:6) = -stable_evec(4:6);

        % Step into manifold
        x_manifold_u_p = xout(i,1:6)' + d * unstable_evec;
        x_manifold_u_n = xout(i,1:6)' - d * unstable_evec;

        % If x-velocity is positive, moon-bound
        % If x-velocity if negative, earth-bound
        if (x_manifold_u_p(4) > 0)
            moon_unstable = x_manifold_u_p;
            earth_unstable = x_manifold_u_n;
        else
            moon_unstable = x_manifold_u_n;
            earth_unstable = x_manifold_u_p;
        end
        
        % Propagate using the event functions
        [moon_unstable_t, moon_unstable_x] = ode113(@(t, state)CR3BP(state, mu), [0, manifold_time], moon_unstable, options);
        [earth_unstable_t, earth_unstable_x] = ode113(@(t, state)CR3BP(state, mu), [0, manifold_time], earth_unstable, options);

        % plot(moon_stable_x(:,1), moon_stable_x(:,2), 'blue')
        % plot(earth_stable_x(:,1), earth_stable_x(:,2), 'red')

        if (abs(moon_unstable_x(end,1) - (1-mu)) < 1e-6 && moon_unstable_x(end,2) > 0)
            moon_unstable_cnt = moon_unstable_cnt + 1;
            moon_bound_stable(:,moon_unstable_cnt) = moon_unstable;
        else
            moon_unstable_cnt = moon_unstable_cnt + 1;
            moon_bound_stable(:,moon_unstable_cnt) = earth_unstable;
        end
    end

    % global count;
    % global poincare_stored;
    % poincare_stored = [];
    % for k = 1:moon_stable_cnt
    %     count = 0;
    %     options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,y) c_eventFn(t, y, mu, n_crossings));
    %     [moon_stable_t, moon_stable_x] = ode113(@(t, state)CR3BP(state, mu), [0, -manifold_time], moon_bound_stable(:,k), options);
    %     plot(moon_stable_x(:,1), moon_stable_x(:,2), 'blue')
    % end
    % hold off
    % legend("Lyapunov Orbit", "L1", "Earth", "Moon")
    % grid on
    % axis equal
    % xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    % ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)    
end

function [value,isterminal,direction] = c_eventFn(t,y,mu,n_crossings)
    global count;
    global poincare_stored;
    global counted_once;
    if n_crossings == 1
        % difference = y(1) - (1-mu);        
        % if abs(difference) < 1e-12
        %     value = 0;
        % else
        %     value = 1;
        % end
        value = y(1) - (1-mu);
        isterminal = 1;
        direction = -1;
        if (abs(value) < 1e-12 && y(4) < 0 && counted_once == false)
            poincare_stored = [poincare_stored; y(2), y(5)];
            counted_once = true;
        end
    else
        if count < n_crossings
            value = y(1) - (1-mu);
            isterminal = 0;
            direction = -1;
            if (abs(value) < 1e-12 && y(4) < 0)
                count = count + 1;
                poincare_stored = [poincare_stored; y(2), y(5)];
                
            end
        elseif count == n_crossings
            value = y(1) - (1-mu); % Want x to be 1-mu
            isterminal = 1; % Halt integration when value is 0
            direction = -1; % When zero is approached from +ve i.e. x_dot > 0
            if (abs(value) < 1e-12 && y(4) < 0)
                poincare_stored = [poincare_stored; y(2), y(5)];
            end
        end
    end
end

function [value,isterminal,direction] = b_eventFn(t,y,mu,n_crossings)
    global count;
    global poincare_stored;
    global counted_once;
    if n_crossings == 1
        % difference = y(1) - (1-mu);        
        % if abs(difference) < 1e-8
        %     value = 0;
        % else
        %     value = 1;
        % end
        value = y(1) - (1-mu);
        isterminal = 1;
        direction = 1;
        if (abs(value) < 1e-12 && y(4) < 0 && counted_once == false)
            poincare_stored = [poincare_stored; y(2), y(5)];    
            counted_once = true;
        end
    else
        if count < n_crossings
            value = y(1) - (1-mu);
            isterminal = 0;
            direction = 1;
            if (abs(value) < 1e-12 && y(4) < 0)
                count = count + 1;
                poincare_stored = [poincare_stored; y(2), y(5)];
                
            end
        elseif count == n_crossings
            value = y(1) - (1-mu); % Want x to be 1-mu
            isterminal = 1; % Halt integration when value is 0
            direction = 1; % When zero is approached from +ve i.e. x_dot > 0
            if (abs(value) < 1e-12 && y(4) < 0)
                poincare_stored = [poincare_stored; y(2), y(5)]; 
            end
        end
    end
end

function [value,isterminal,direction] = eventFn(t,y, mu)
    value = [1-mu-y(1), y(1)-(-mu)];
    isterminal = [1, 1]; % Halt integration when value is 0
    direction = [0, 0]; % When zero is approached from either side
end

function out = isnotin(array, val)
    out = true;
    for el = 1:length(array)
        if val == array(el)
            out = false;
        end
    end
end


