clear; clc; close all;

%% ASEN 5090 - HW 1
% Jash Bhalavat
% 08/24/2026

% Constants

% Speed of light in a vacuum
c = 2.99792458e8; % m/s, ME Table 4.1

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
