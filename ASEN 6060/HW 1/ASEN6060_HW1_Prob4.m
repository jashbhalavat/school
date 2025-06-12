clear; clc; close all;

% ASEN 6060 - HW 1, Problem 4
% Jash Bhalavat
% 01/28/2025

%% Constants
G = 6.67408 * 10^-11; % m3/(kgs2)
G = G / (10^9); % km3/(kgs2)
% Earth
mu_earth = 398600.435507; % km3/s2
mass_earth = mu_earth / G; % kg
% Moon
mu_moon = 4902.800118; % km3/s2
a_moon = 384400; % km
mass_moon = mu_moon / G; % kg
% Earth-Moon system
mass_ratio_em = mass_moon / (mass_earth + mass_moon);
m_star_em = mass_earth + mass_moon;
mu = mass_ratio_em;

%% Zero Velocity Surfaces

mu = mass_ratio_em;

n = 1001;
x = linspace(-2, 2, n);
y = linspace(-2, 2, n);

% Given Jacobi constant
c_case = 4;

switch c_case
    case 1
        c_given = 3.189;
    case 2
        c_given = 3.173;
    case 3
        c_given = 3.013;
    case 4
        c_given = 2.995;
end

zvc = full_zvc(c_given, mu, x, y);
plot_zvc(zvc, x, y, c_given, mu)

%% Part d

c_test = linspace(2,3.2,1000);

% Y index for zero
y_zero_ind = find(y == 0); 

% Find -mu and 1 - mu X indices
[del, x_p1_ind] = min(abs(x - (-mu)));
[del, x_p2_ind] = min(abs(x - (1-mu)));

l1_calc = false;
l2_calc = false;
l3_calc = false;
l45_calc = false;

for i = 1:length(c_test)
    zvc = full_zvc(c_test(i), mu, x, y);

    % Check for C(L4/5)
    % Check if zvc has any dynamically excluded motion
    zvc_neg_1 = find(zvc == -1);
    if and(~isempty(zvc_neg_1), ~l45_calc)
        c_l45 = c_test(i-1);
        zvc_l45 = full_zvc(c_l45, mu, x, y);
        zvc_l45_minus_1 = full_zvc(c_test(i), mu, x, y);
        l45_calc = true;
    end

    % Check for C(L3)
    % L3 is on y = 0
    % Check where dynamically excluded motion first appears in that row
    zvc_zero_y = zvc(:, y_zero_ind);
    zvc_zero_y_neg_1 = find(zvc_zero_y == -1);
    if and(~isempty(zvc_zero_y_neg_1), ~l3_calc)
        c_l3 = c_test(i);
        zvc_l3 = full_zvc(c_l3, mu, x, y);
        l3_calc = true;
    end

    % Check for C(l2)
    % L2 is located beyond P2 on y=0
    zvc_beyond_p2 = zvc(x_p2_ind:end, y_zero_ind);
    zvc_beyond_p2_neg_1 = find(zvc_beyond_p2 == -1);
    if and(~isempty(zvc_beyond_p2_neg_1), ~l2_calc)
        c_l2 = c_test(i);
        zvc_l2 = full_zvc(c_l2, mu, x, y);
        l2_calc = true;
    end

    % Check for C(l1)
    % L1 is located between P1 and P2 on y=0
    zvc_between_p1_p2 = zvc(x_p1_ind:x_p2_ind, y_zero_ind);
    zvc_between_p1_p2_neg_1 = find(zvc_between_p1_p2 == -1);
    if and(~isempty(zvc_between_p1_p2_neg_1), ~l1_calc)
        c_l1 = c_test(i);
        zvc_l1 = full_zvc(c_l1, mu, x, y);
        l1_calc = true;
    end

end

plot_zvc(zvc_l1, x, y, c_l1, mu)
plot_zvc(zvc_l2, x, y, c_l2, mu)
plot_zvc(zvc_l3, x, y, c_l3, mu)
plot_zvc(zvc_l45, x, y, c_l45, mu)
plot_zvc(zvc_l45_minus_1, x, y, c_l45, mu)

function zvc = full_zvc(c_given, mu, x, y)
    zvc = ones([length(x), length(y)]);
    
    for i = 1:length(x)
        for j = 1:length(y)
            c_calc = u_star_times_2(x(i), y(j), mu);
            if c_calc < c_given
                zvc(i, j) = -1;
            end
        end
    end
end

function plot_zvc(zvc, x, y, c, mu)
    figure()
    contourf(x, y, zvc')
    map = [220/255, 220/255, 220/255
        1, 1, 1];
    colormap(map)
    hold on
    scatter(-mu, 0, 200, 'filled')
    scatter(1-mu, 0, 50, 'filled')
    hold off
    legend("", "Earth", "Moon")
    xlabel("x [Non-Dimensional]")
    ylabel("y [Non-Dimensional]")
    title("Zero velocity curve for Jacobi Constant C = " + num2str(c))
end

function out = u_star_times_2(x, y, mu)
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2);
    out = (x^2 + y^2) + 2*(1 - mu)/r1 + 2*mu/r2;
end


