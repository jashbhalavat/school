clear; clc; close all;

% ASEN 6060 - HW 2, Problem 1
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

%% Part a, b

mu = mass_ratio_em;

% Earth Moon system equilibrium points
[em_eq_pts, em_eq_validity] = all_eq_points(mu);

% Jacobi constant at each eq point
for i = 1:5
    x = em_eq_pts(i, 1);
    y = em_eq_pts(i, 2);
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2);
    c_at_eq(i) = (x^2 + y^2) + 2*(1 - mu)/r1 + 2*mu/r2;
end

%% Part c

mu_range = linspace(1e-7, 0.5, 100);

for i = 1:length(mu_range)
    [eq_pts(:,:,i), eq_validity(:,:,i)] = all_eq_points(mu_range(i));
end

figure()
subplot(2,2,1)
plot(mu_range, squeeze(eq_pts(1,1,:)), 'LineWidth',2)
xlabel("Mass Ratio [mu]")
ylabel("L1 x position")

subplot(2,2,2)
plot(mu_range, squeeze(eq_pts(2,1,:)), 'LineWidth',2)
xlabel("Mass Ratio [mu]")
ylabel("L2 x position")

subplot(2,2,3)
plot(mu_range, squeeze(eq_pts(3,1,:)), 'LineWidth',2)
xlabel("Mass Ratio [mu]")
ylabel("L3 x position")

subplot(2,2,4)
plot(mu_range, squeeze(eq_pts(4,1,:)), 'LineWidth',2)
xlabel("Mass Ratio [mu]")
ylabel("L4/5 x position")

sgtitle("Equilibrium Points x position as a function of mass ratio")


