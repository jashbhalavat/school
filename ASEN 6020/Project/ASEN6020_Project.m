clear; clc; close all;

% ASEN 6020 - Project
% Spring 2025
% Jash Bhalavat

global angle_i delay_time_i duration_i

%% Angle

options = optimoptions('fmincon', 'ConstraintTolerance',1e-12);
angle_i = [];
angle_star = fmincon(@cost_fn, 0.5*pi/180, [], [], [], [], -90*pi/180, 90*pi/180, [], options);

one_over_delta_i_star = cost_fn(angle_star);
delta_i_star = 1/one_over_delta_i_star;

figure()
plot(1:length(angle_i), angle_i*180/pi, 'LineWidth', 2)
grid on
xlabel("Iterations")
ylabel("Optimal Thrust Angle [°]")
title("Optimal Thrust Angle when thrusting for 1 full period")

%% Delay Time

r_0 = 6928; % km
mu = 3.986004415e5;
P = 2*pi*sqrt(r_0^3/mu);

delay_time_i = [];
delay_time_star = fmincon(@cost_fn_delay, 100, [], [], [], [], 0, P/2-1, [], options);

one_over_di_star = cost_fn(delay_time_star);
di_star = 1/one_over_di_star;

figure()
plot(1:length(delay_time_i), delay_time_i, 'LineWidth', 2)
grid on
xlabel("Iterations")
ylabel("Delay Time [s]")
title("Optimal Delay Time")

%% Fire Duration


duration_i = [];
duration_star = fmincon(@cost_fn_duration, 250, [], [], [], [], 0, P/2-1, [], options);

one_over_di_star = cost_fn(duration_star);
di_star = 1/one_over_di_star;

figure()
plot(1:length(duration_i), duration_i, 'LineWidth', 2)
grid on
xlabel("Iterations")
ylabel("Fire Duration [s]")
title("Optimal Fire Duration")

%% One Period Plots

thrust = 0;
I_sp = 1000;
g_0 = 9.8;
mass_0 = 1000; % kg
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

r_0 = 6928; % km
v_0 = sqrt(mu/r_0);

coe_0 = [r_0, 0, 98.8*pi/180, 0, 0, 0];
rv_0 = coe2rv(coe_0, mu);

P = 2*pi*sqrt(r_0^3/mu);

[tout, xout] = ode45(@(t, state)TBP(state, mu, thrust, I_sp, g_0, 0), [0, P], [rv_0; mass_0], options);

figure()
subplot(2,1,1)
plot(tout, xout(:,1:3), 'LineWidth',2)
legend("x", "y", "z")
grid on
title("No Thrust Orbit")
xlabel("Time [sec]")
ylabel("Inertial Position [km]")

thrust = 100 * 10^-3;
angle_star = 0.800307432766832;
delay_time_star = 1.430447642010037e+03;

[tout_star, xout_star] = ode45(@(t, state)TBP_duration(t, state, mu, thrust, I_sp, g_0, angle_star, delay_time_star, duration_star), [0, P], [rv_0; mass_0], options);

subplot(2,1,2)
plot(tout_star, xout_star(:,1:3), 'LineWidth',2)
grid on
title("Orbit w/Thrust")
xregion(delay_time_star, duration_star, FaceColor='blue')
legend("x", "y", "z", "Thrusting Region")
xlabel("Time [sec]")
ylabel("Inertial Position [km]")

sgtitle("Spacecraft Inertial Position States")





function one_over_delta_i = cost_fn_duration(duration)
    global duration_i
    duration_i = [duration_i, duration];
    mu = 3.986004415e5;
    thrust = 100 * 10^-3;
    I_sp = 1000;
    g_0 = 9.8;
    mass_0 = 1000; % kg
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    r_0 = 6928; % km
    v_0 = sqrt(mu/r_0);
    
    coe_0 = [r_0, 0, 98.8*pi/180, 0, 0, 0];
    rv_0 = coe2rv(coe_0, mu);
    
    P = 2*pi*sqrt(r_0^3/mu);

    angle_star = 0.800307432766832;
    delay_time_star = 1.430447642010037e+03;
    
    [tout, xout] = ode45(@(t, state)TBP_duration(t, state, mu, thrust, I_sp, g_0, angle_star, delay_time_star, duration), [0, P], [rv_0; mass_0], options);
    
    state_f = xout(end,:);
    coe_f = rv2coe(state_f(1:3), state_f(4:6), mu);
    delta_i = (coe_f(3) - coe_0(3)) * 180/pi;
    one_over_delta_i = abs(1/delta_i);
end


function one_over_delta_i = cost_fn_delay(delay_time)
    global delay_time_i
    delay_time_i = [delay_time_i, delay_time];
    mu = 3.986004415e5;
    thrust = 100 * 10^-3;
    I_sp = 1000;
    g_0 = 9.8;
    mass_0 = 1000; % kg
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    r_0 = 6928; % km
    v_0 = sqrt(mu/r_0);
    
    coe_0 = [r_0, 0, 98.8*pi/180, 0, 0, 0];
    rv_0 = coe2rv(coe_0, mu);
    
    P = 2*pi*sqrt(r_0^3/mu);

    angle_star = 0.800307432766832;
    
    [tout, xout] = ode45(@(t, state)TBP_delay_time(t, state, mu, thrust, I_sp, g_0, angle_star, delay_time, P), [0, P], [rv_0; mass_0], options);
    
    state_f = xout(end,:);
    coe_f = rv2coe(state_f(1:3), state_f(4:6), mu);
    delta_i = (coe_f(3) - coe_0(3)) * 180/pi;
    one_over_delta_i = abs(1/delta_i);
end

function one_over_delta_i = cost_fn(angle)
    global angle_i
    angle_i = [angle_i, angle];
    fprintf("Angle - %d\n", angle*180/pi)
    mu = 3.986004415e5;
    thrust = 100 * 10^-3;
    I_sp = 1000;
    g_0 = 9.8;
    mass_0 = 1000; % kg
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    r_0 = 6928; % km
    v_0 = sqrt(mu/r_0);
    
    coe_0 = [r_0, 0, 98.8*pi/180, 0, 0, 0];
    rv_0 = coe2rv(coe_0, mu);
    
    P = 2*pi*sqrt(r_0^3/mu);
    
    [tout, xout] = ode45(@(t, state)TBP(state, mu, thrust, I_sp, g_0, angle), [0, P], [rv_0; mass_0], options);
    
    state_f = xout(end,:);
    coe_f = rv2coe(state_f(1:3), state_f(4:6), mu);
    delta_i = (coe_f(3) - coe_0(3)) * 180/pi;
    one_over_delta_i = abs(1/delta_i);
end
% dv = sqrt((state_f(4) - rv_0(4))^2 + (state_f(5) - rv_0(5))^2 + (state_f(6) - rv_0(6))^2);

