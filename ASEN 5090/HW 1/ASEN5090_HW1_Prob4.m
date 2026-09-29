clear; clc; close all;

%% ASEN 5090 - HW 1
% Jash Bhalavat
% 08/25/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

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