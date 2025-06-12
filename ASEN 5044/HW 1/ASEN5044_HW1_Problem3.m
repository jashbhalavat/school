clear; clc; close all;

syms x

p0 = 20;
Ix = 500;
Iy = 750;
Iz = 1000;

A = [0 0 0; 0 0 x*p0*(Ix-Iz)/Iy; 0 x*p0*(Iy-Ix)/Iz 0];

stm = simplify(expm(A));
stmf(x) = stm;
stmf(0.1)

% Part c
time = linspace(0, 5, 1000);

x0 = [0 0.1 0]';

for i = 1:length(time)
    system(:, i) = stmf(time(i)) * x0;
end

figure()
plot(time, system(1, :))
hold on
plot(time, system(2, :))
plot(time, system(3, :))
hold off
legend("Rolling", "Pitching", "Yawing")
xlabel("Time [s]")
ylabel("Perturbation [rad]")
title("State Time History")
