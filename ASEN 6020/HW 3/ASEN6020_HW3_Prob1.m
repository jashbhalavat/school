clear; clc; close all;

% ASEN 6020 - HW 3 - Problem 1b
% Spring 2025
% Jash Bhalavat

%% Part b

global_x0 = [0;1;-0.5];
global_root = fsolve(@(x)fx(x), global_x0);
global_root_y2 = 2 + sin(pi*global_root(3));
d_global = sqrt((global_root(1)-global_root(3))^2 + (global_root(2)-global_root_y2)^2);

% first_x0 = [-0.5;1];
first_x0 = [0;1;-2.5];
first_root = fsolve(@(x)fx(x), first_x0);
first_root_y2 = 2 + sin(pi*first_root(3));
d_first = sqrt((first_root(1)-first_root(3))^2 + (first_root(2)-first_root_y2)^2);

% second_x0 = [1.5;1];
second_x0 = [0;1;1.5];
second_root = fsolve(@(x)fx(x), second_x0);
second_root_y2 = 2 + sin(pi*second_root(3));
d_second = sqrt((second_root(1)-second_root(3))^2 + (second_root(2)-second_root_y2)^2);

% third_x0 = [-2.5;1];
third_x0 = [0;1;3.5];
third_root = fsolve(@(x)fx(x), third_x0);
third_root_y2 = 2 + sin(pi*third_root(3));
d_third = sqrt((third_root(1)-third_root(3))^2 + (third_root(2)-third_root_y2)^2);

x1 = linspace(-10, 10, 1000);
for i = 1:length(x1)
    y1(i) = sqrt(1 - x1(i)^2);
    y2(i) = 2 + sin(pi*x1(i));
end

figure()
plot(x1, y1)
hold on
plot(x1, y2)
scatter(global_root(1), global_root(2), 'filled', 'magenta')
scatter(global_root(3), global_root_y2, 'filled', 'cyan')
% scatter(first_root(1), first_root(2), 'filled', 'black')
scatter(first_root(3), first_root_y2, 'filled', 'black')
% scatter(second_root(1), second_root(2), 'filled', 'blue')
scatter(second_root(3), second_root_y2, 'filled', 'blue')
% scatter(third_root(1), third_root(2), 'filled', 'red')
scatter(third_root(3), third_root_y2, 'filled', 'red')
legend("g1", "g2", "x1, y1 for all cases","Global (x2, y2)", "Case 1 (x2, y2)", "Case 2 (x2, y2)", "Case 3 (x2, y2")
grid on
title("Terminal Manifolds with Global and Local Minima")
xlabel("x")
ylabel("y")

function out = fx(x)
    out = x(1) + pi * x(2) * cos(pi*x(3));
end