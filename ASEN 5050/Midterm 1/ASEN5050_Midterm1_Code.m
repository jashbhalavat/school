clear; clc; close all;

% ASEN 5050, Midterm 1
% Fall 20204
% Author - Jash Bhalavat

GM_mars = 4.305e4;
r_eq_mars = 3397.2;

%% Problem 1

h = 2.4799e4;

% At t1, phi_fpa < 0
r1 = 15000;
v1 = 1.8111;

% At t2, v_r2 < 0
r2 = 19000;

% Part a
eps = v1^2/2 - GM_mars/r1;
a = -GM_mars/(2*eps);

% Part b
e = sqrt((2*h^2 * eps)/(GM_mars^2) + 1);

% Part c
E1 = acos(1/e*(1 - r1/a)) * -1;
theta_star_1 = acos(((h^2/GM_mars)/r1 - 1)/e) * -1;

% Part d
E2 = acos(1/e*(1 - r2/a)) * -1;
theta_star_2 = acos(((h^2/GM_mars)/r2 - 1)/e) * -1;

% Part e
P = 2*pi * sqrt(a^3/GM_mars);
n = sqrt(GM_mars/a^3);

t1_to_tp = 1/n * abs(E1 - e*sin(E1));
t2_to_tp = 1/n * abs(E2 - e*sin(E2));
tp_to_t2 = P - t2_to_tp;

t1_to_t2 = t1_to_tp + tp_to_t2;

% Part f
rp = a * (1 - e);
ra = a * (1 + e);
b = a * sqrt(1 - e^2);
p = a * (1 - e^2);


%% Problem 2

% At time t3
r3_relay = [-7.6650e3, 6.5468e3, -4.5740e2];
v3_relay = [1.6334, 0.1226, -1.9455];

% Part a
coe_t3_relay = rv2coe(r3_relay, v3_relay, GM_mars);
theta_star_3_relay = coe_t3_relay(6);
e_relay = coe_t3_relay(2);
a_relay = coe_t3_relay(1);

% Part b
% Time t4 is t3 + 2 hours
P_relay = 2*pi * sqrt(a_relay^3/GM_mars);
n_relay = sqrt(GM_mars/a_relay^3);

E3_relay = 2 * atan(sqrt((1-e_relay)/(1+e_relay)) * tan(theta_star_3_relay/2));
t3_to_tp = 1/n_relay * (E3_relay - e_relay*sin(E3_relay));

% T3 is before peiapsis
% T4 is time past periapsis
t4 = t3_to_tp + 2*3600;

E4_relay = kepler_solver_eclipse(t4, a_relay, e_relay, GM_mars);
theta_star_4_relay = atan(sqrt((1 + e_relay)/(1 - e_relay)) * tan(E4_relay/2)) * 2;

r4v4_relay = fg(r3_relay, v3_relay, theta_star_4_relay, GM_mars);


%% Problem 3

r_cubesat_hat_descending_node = [-0.64279, -0.76604, 0];
r_cubesat_hat_apoapsis = [-0.02970, -0.97508, -0.21985];

n_hat_cubesat = -1 * r_cubesat_hat_descending_node;
e_hat_cubesat = -1 * r_cubesat_hat_apoapsis;

raan = sign(n_hat_cubesat(2)) * abs(acos(n_hat_cubesat(1)));
aop = sign(e_hat_cubesat(3)) * abs(acos(dot(n_hat_cubesat , e_hat_cubesat)));









