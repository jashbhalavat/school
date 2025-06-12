clear; clc; close all;

% ASEN 5010 - HW 4, Problem 2
% Spring 2025
% Jash Bhalavat

I = [125, 100, 75];
L = [0, 0, 0]';

state1 = [1, 0, 0, 0, 0, 0]';
state2 = [0, 1, 0, 0, 0, 0]';
state3 = [0, 0, 1, 0, 0, 0]';
state4 = [1, 0.1, 0, 0, 0, 0]';
state5 = [1, 0, 0.1, 0, 0, 0]';
state6 = [0.1, 1, 0, 0, 0, 0]';
state7 = [0, 1, 0.1, 0, 0, 0]';
state8 = [0, 0.1, 1, 0, 0, 0]';
state9 = [0.1, 0, 1, 0, 0, 0]';

tspan = [0, 100];

[tout1, state_out1] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state1);
[tout2, state_out2] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state2);
[tout3, state_out3] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state3);
[tout4, state_out4] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state4);
[tout5, state_out5] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state5);
[tout6, state_out6] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state6);
[tout7, state_out7] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state7);
[tout8, state_out8] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state8);
[tout9, state_out9] = ode45(@(t, y) diff_eq(t, y, I, L), tspan, state9);

figure(1)
subplot(3,2,1)
plot(tout1, state_out1(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state1(1)), ', ', num2str(state1(2)), ', ', num2str(state1(3)), '] rad/s'])
subplot(3,2,3)
plot(tout2, state_out2(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state2(1)), ', ', num2str(state2(2)), ', ', num2str(state2(3)), '] rad/s'])
subplot(3,2,5)
plot(tout3, state_out3(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state3(1)), ', ', num2str(state3(2)), ', ', num2str(state3(3)), '] rad/s'])

subplot(3,2,2)
plot(tout1, state_out1(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
title(['\omega_0 = [', num2str(state1(1)), ', ', num2str(state1(2)), ', ', num2str(state1(3)), '] rad/s'])
subplot(3,2,4)
plot(tout2, state_out2(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
title(['\omega_0 = [', num2str(state2(1)), ', ', num2str(state2(2)), ', ', num2str(state2(3)), '] rad/s'])
subplot(3,2,6)
plot(tout3, state_out3(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
title(['\omega_0 = [', num2str(state3(1)), ', ', num2str(state3(2)), ', ', num2str(state3(3)), '] rad/s'])
sgtitle("Initial body rates 1 rad/s about each principal axis")

figure(2)
subplot(2,2,1)
plot(tout4, state_out4(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state4(1)), ', ', num2str(state4(2)), ', ', num2str(state4(3)), '] rad/s'])
subplot(2,2,3)
plot(tout5, state_out5(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state5(1)), ', ', num2str(state5(2)), ', ', num2str(state5(3)), '] rad/s'])

subplot(2,2,2)
plot(tout4, state_out4(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
ylim([-1/2 1/2])
title(['\omega_0 = [', num2str(state4(1)), ', ', num2str(state4(2)), ', ', num2str(state4(3)), '] rad/s'])
subplot(2,2,4)
plot(tout5, state_out5(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
ylim([-1/2 1/2])
title(['\omega_0 = [', num2str(state5(1)), ', ', num2str(state5(2)), ', ', num2str(state5(3)), '] rad/s'])
sgtitle("Initial body rates 1 rad/s about each principal axis and 0.1 rad/s about another axis")

figure(3)
subplot(2,2,1)
plot(tout6, state_out6(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state6(1)), ', ', num2str(state6(2)), ', ', num2str(state6(3)), '] rad/s'])
subplot(2,2,3)
plot(tout7, state_out7(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state7(1)), ', ', num2str(state7(2)), ', ', num2str(state7(3)), '] rad/s'])

subplot(2,2,2)
plot(tout6, state_out6(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
title(['\omega_0 = [', num2str(state6(1)), ', ', num2str(state6(2)), ', ', num2str(state6(3)), '] rad/s'])
subplot(2,2,4)
plot(tout7, state_out7(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
title(['\omega_0 = [', num2str(state7(1)), ', ', num2str(state7(2)), ', ', num2str(state7(3)), '] rad/s'])
sgtitle("Initial body rates 1 rad/s about each principal axis and 0.1 rad/s about another axis")

figure(4)
subplot(2,2,1)
plot(tout8, state_out8(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state8(1)), ', ', num2str(state8(2)), ', ', num2str(state8(3)), '] rad/s'])
subplot(2,2,3)
plot(tout9, state_out9(:,1:3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("\omega [rad/s]")
legend("b1", "b2", "b3")
ylim([-2 2])
title(['\omega_0 = [', num2str(state9(1)), ', ', num2str(state9(2)), ', ', num2str(state9(3)), '] rad/s'])

subplot(2,2,2)
plot(tout8, state_out8(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
title(['\omega_0 = [', num2str(state8(1)), ', ', num2str(state8(2)), ', ', num2str(state8(3)), '] rad/s'])
subplot(2,2,4)
plot(tout9, state_out9(:,4:6), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Euler Angle [rad]")
legend("\psi", "\theta", "\phi")
title(['\omega_0 = [', num2str(state9(1)), ', ', num2str(state9(2)), ', ', num2str(state9(3)), '] rad/s'])
sgtitle("Initial body rates 1 rad/s about each principal axis and 0.1 rad/s about another axis")


% function ea_dot = ea_diff_eq(t, ea, omega)
%     psi = ea(1);
%     theta = ea(2);
%     phi = ea(3);
%     ea_dot = 1/cos(theta) * [0, sin(phi), cos(phi);
%                             0, cos(phi)*cos(theta), -sin(phi)*cos(theta);
%                             cos(theta), sin(phi)*sin(theta), cos(phi)*sin(theta)] * omega;
% end


function ydot = diff_eq(t, y, I, L)
    % State - [omega; ea]
    
    ydot(1,1) = -(I(3) - I(2))/I(1) * y(2) * y(3) + L(1);
    ydot(2,1) = -(I(1) - I(3))/I(2) * y(3) * y(1) + L(2);
    ydot(3,1) = -(I(2) - I(1))/I(3) * y(2) * y(1) + L(3);

    psi = y(4);
    theta = y(5);
    phi = y(6);
    ydot(4:6,1) = 1/cos(theta) * [0, sin(phi), cos(phi);
                            0, cos(phi)*cos(theta), -sin(phi)*cos(theta);
                            cos(theta), sin(phi)*sin(theta), cos(phi)*sin(theta)] * y(1:3);
end
