clear; clc; close all;

%% ASEN 5090 - HW 1
% Jash Bhalavat
% 08/25/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

%% Problem 1a

rho_PL1 = 550; % m
rho_PL2 = 500; % m

delta_r = (rho_PL1 - rho_PL2) / (2*c); % sec
R_r_PL1 = rho_PL1 - c*delta_r; % m
R_r_PL2 = 1000 - R_r_PL1; % m

n = 500;

%% Problem 2

v0 = 50; % m/s
h0 = 100; % m
x0 = 250; % m

xA = linspace(-x0, x0, n);

h0_err = 10/100;
h0_neg = h0 - h0_err;
h0_pos = h0 + h0_err;

for i = 1:n
    rho(i) = sqrt(h0^2 + xA(i)^2); % m
    rho_neg(i) = sqrt(h0_neg^2 + xA(i)^2); % m
    rho_pos(i) = sqrt(h0_pos^2 + xA(i)^2); % m
    rho_dot(i) = v0*xA(i) / rho(i); % m/s
    rho_dot_neg(i) = v0*xA(i) / rho_neg(i); % m/s
    rho_dot_pos(i) = v0*xA(i) / rho_pos(i); % m/s
    % z(i) = asin(xA(i) / rho(i)); % 
    z(i) = atan2(xA(i), h0);
    z_neg(i) = atan2(xA(i), h0_neg);
    z_pos(i) = atan2(xA(i), h0_pos);
end

% Plot
figure(1)
subplot(3,1,1)
plot(xA, rho, 'LineWidth', 2)
hold on
plot(xA, rho_neg,  '--', 'Color', 'red', 'LineWidth', 2)
plot(xA, rho_pos, '--', 'Color', 'red', 'LineWidth', 2)
hold off
ylabel('$\rho$ [m]', 'Interpreter', 'Latex', 'FontSize', 15)
grid on

subplot(3,1,2)
plot(xA, rho_dot, 'LineWidth', 2)
hold on
plot(xA, rho_dot_neg,  '--', 'Color', 'red', 'LineWidth', 2)
plot(xA, rho_dot_pos, '--', 'Color', 'red', 'LineWidth', 2)
hold off
ylabel('$\dot{\rho}$ [m/s]', 'Interpreter', 'Latex', 'FontSize', 15)
grid on

subplot(3,1,3)
plot(xA, z, 'LineWidth', 2)
hold on
plot(xA, z_neg,  '--', 'Color', 'red', 'LineWidth', 2)
plot(xA, z_pos, '--', 'Color', 'red', 'LineWidth', 2)
hold off
ylabel('$z$ [rad]', 'Interpreter', 'Latex', 'FontSize', 15)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15)
grid on

sgtitle("Range, Range Rate, and Zenith Angle Measurements")


%% Problem 3a

% For PRN 19, S1 = 3, S2 = 6. Source = https://archive.gps.gov/technical/icwg/IS-GPS-200D.pdf
S1 = 3;
S2 = 6;

ca_code_prn19_1023 = ca_code(1023, S1, S2);

prn19_first16_hex = dec2hex(bin2dec(num2str(ca_code_prn19_1023(1:16))));
prn19_last16_hex = dec2hex(bin2dec(num2str(ca_code_prn19_1023(end-15:end))));

figure(1)
subplot(2,1,1)
stairs(ca_code_prn19_1023(1:16), 'linewidth', 2)
title("First 16 chips")
grid on

subplot(2,1,2)
stairs(ca_code_prn19_1023(end-15:end), 'linewidth', 2)
title("Last 16 chips")
grid on

sgtitle("First and last 16 chips of PRN 19 C/A code")

%% Problem 3b

ca_code_prn19_2046 = ca_code(2046, S1, S2);

figure(2)
stairs(ca_code_prn19_2046(1:1023), 'linewidth', 1)
hold on
stairs(ca_code_prn19_2046(1024:2046), 'linewidth', 1)
hold off
legend("Epoch 1-1023", "Epoch 1024-2046")
grid on

title("PRN 19 C/A code for Epoch 1-1023 and 1024-2046")

%% Problem 3c

% For PRN 25, S1 = 5, S2 = 7. Source = https://archive.gps.gov/technical/icwg/IS-GPS-200D.pdf
S1 = 5;
S2 = 7;

ca_code_prn25_1023 = ca_code(1023, S1, S2);

prn25_first16_hex = dec2hex(bin2dec(num2str(ca_code_prn25_1023(1:16))));
prn25_last16_hex = dec2hex(bin2dec(num2str(ca_code_prn25_1023(end-15:end))));

figure(3)
subplot(2,1,1)
stairs(ca_code_prn25_1023(1:16), 'linewidth', 2)
title("First 16 chips")
grid on

subplot(2,1,2)
stairs(ca_code_prn25_1023(end-15:end), 'linewidth', 2)
title("Last 16 chips")
grid on

sgtitle("First and last 16 chips of PRN 25 C/A code")

%% Problem 3d

% For PRN 5, S1 = 1, S2 = 9. Source = https://archive.gps.gov/technical/icwg/IS-GPS-200D.pdf
S1 = 1;
S2 = 9;

ca_code_prn5_1023 = ca_code(1023, S1, S2);

prn5_first16_hex = dec2hex(bin2dec(num2str(ca_code_prn5_1023(1:16))));
prn5_last16_hex = dec2hex(bin2dec(num2str(ca_code_prn5_1023(end-15:end))));

figure(4)
subplot(2,1,1)
stairs(ca_code_prn5_1023(1:16), 'linewidth', 2)
title("First 16 chips")
grid on

subplot(2,1,2)
stairs(ca_code_prn5_1023(end-15:end), 'linewidth', 2)
title("Last 16 chips")
grid on

sgtitle("First and last 16 chips of PRN 5 C/A code")

%% Problem 4a

% For PRN 19, S1 = 3, S2 = 6. Source = https://archive.gps.gov/technical/icwg/IS-GPS-200D.pdf
S1 = 3;
S2 = 6;

ca_code_prn19_1023 = ca_code(1023, S1, S2);
ca_problem2_prn19_1023 = (ca_code_prn19_1023==0)*(1) + (ca_code_prn19_1023==1)*(-1);

n = linspace(-1023, 1023, 2047);

for i = 1:length(n)
    R19(i) = auto_correlation_norm(n(i), ca_problem2_prn19_1023);
end

figure(1)
plot(n, R19, 'LineWidth', 1.5)
xlim([-1050, 1050])
xlabel("Shift (n)")
ylabel("$R^{19}$", 'Interpreter', 'latex')
title("PRN19 Normalized Auto-Correlation as a function of Shift")
grid on
ylim([-0.2, 1.2])


%% Problem 4b

ca_problem2_prn19_1023_shifted_by_200 = circshift(ca_problem2_prn19_1023, [0, 200]);

for i = 1:length(n)
    R19_19200(i) = cross_correlation_norm(n(i), ca_problem2_prn19_1023_shifted_by_200, ca_problem2_prn19_1023);
end

figure(2)
plot(n, R19_19200, 'LineWidth', 1.5)
xlim([-1050, 1050])
xlabel("Shift (n)")
ylabel("$R^{19,19_{200}}$", 'Interpreter', 'latex')
title("Normalized Cross-Correlation of PRN19 and PRN19 shifted by 200 chips")
grid on
ylim([-0.2, 1.2])

%% Problem 4c

% For PRN 25, S1 = 5, S2 = 7. Source = https://archive.gps.gov/technical/icwg/IS-GPS-200D.pdf
S1 = 5;
S2 = 7;

ca_code_prn25_1023 = ca_code(1023, S1, S2);
ca_problem2_prn25_1023 = (ca_code_prn25_1023==0)*(1) + (ca_code_prn25_1023==1)*(-1);

for i = 1:length(n)
    R19_25(i) = cross_correlation_norm(n(i), ca_problem2_prn25_1023, ca_problem2_prn19_1023);
end

figure(3)
plot(n, R19_25, 'LineWidth', 1.5)
xlim([-1050, 1050])
xlabel("Shift (n)")
ylabel("$R^{19,25}$", 'Interpreter', 'latex')
title("Normalized Cross-Correlation of PRN19 and PRN25")
grid on
ylim([-0.2, 1.2])

%% Problem 4d

% For PRN 5, S1 = 1, S2 = 9. Source = https://archive.gps.gov/technical/icwg/IS-GPS-200D.pdf
S1 = 1;
S2 = 9;

ca_code_prn5_1023 = ca_code(1023, S1, S2);
ca_problem2_prn5_1023 = (ca_code_prn5_1023==0)*(1) + (ca_code_prn5_1023==1)*(-1);

for i = 1:length(n)
    R19_5(i) = cross_correlation_norm(n(i), ca_problem2_prn5_1023, ca_problem2_prn19_1023);
end

figure(4)
plot(n, R19_5, 'LineWidth', 1.5)
xlim([-1050, 1050])
xlabel("Shift (n)")
ylabel("$R^{19,5}$", 'Interpreter', 'latex')
title("Normalized Cross-Correlation of PRN19 and PRN5")
grid on
ylim([-0.2, 1.2])

%% Problem 4e

x1 = circshift(ca_problem2_prn19_1023, 350);
x2 = circshift(ca_problem2_prn25_1023, 905);
x3 = circshift(ca_problem2_prn5_1023, 75);
x123 = x1 + x2 + x3;

for i = 1:length(n)
    Rx123_19(i) = cross_correlation_norm(n(i), x123, ca_problem2_prn19_1023);
end

figure(5)
plot(n, Rx123_19, 'LineWidth', 1.5)
xlim([-1050, 1050])
xlabel("Shift (n)")
ylabel("$R^{x123,19}$", 'Interpreter', 'latex')
title("Normalized Cross-Correlation of x1+x2+x3 and PRN19")
grid on
ylim([-0.2, 1.2])

%% Problem 4f

noise = 4*randn(1,1023);

figure(6)
sgtitle("x1, x2, x3 and noise")

subplot(4, 1, 1)
stairs(x1)
ylabel("x1")

subplot(4, 1, 2)
stairs(x2)
ylabel("x2")

subplot(4, 1, 3)
stairs(x3)
ylabel("x3")

subplot(4, 1, 4)
plot(noise)
ylabel("Noise")
xlabel("Chips")

%% Problem 4g

x123n = x123 + noise;

for i = 1:length(n)
    Rx123n_19(i) = cross_correlation_norm(n(i), x123n, ca_problem2_prn19_1023);
end

figure(7)
plot(n, Rx123n_19, 'LineWidth', 1.5)
xlim([-1050, 1050])
xlabel("Shift (n)")
ylabel("$R^{x123n,19}$", 'Interpreter', 'latex')
title("Normalized Cross-Correlation of x1+x2+x3+noise and PRN19")
grid on
ylim([-0.2, 1.2])

%% Functions

function out = auto_correlation_norm(n, x_k)
    % Normalized auto correlation function
    % n - shift
    % x_k - PRN associated with satellite k

    out = sum(x_k .* circshift(x_k, n))/1023;
    
end

function out = cross_correlation_norm(n, x_k, x_l)
    % Normalized cross correlation function
    % n - shift
    % x_k - PRN associated with satellite k
    % x_l - PRN associated with satellite l

    out = sum(x_k .* circshift(x_l, [0, n]))/1023;
    
end