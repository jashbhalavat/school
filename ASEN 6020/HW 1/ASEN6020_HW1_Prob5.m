close all; clear; clc;

% ASEN 6020, HW 1, Problem 5
% Jash Bhalavat
% Spring 2025

delta_i_range = linspace(0.001, pi, 10);
l_range = 1:100;

x0 = 0.5;

tol = 1e-10;

% Set display to zero and thresholds to tol
options = optimoptions('fsolve','Display','none','FunctionTolerance',tol,'OptimalityTolerance',tol,'StepTolerance',tol);

dl_maneuver(0.25, 0.05, 10)

for i = 1:length(delta_i_range)
    for j = 1:length(l_range)
        delta_i = delta_i_range(i);
        l = l_range(j);

        % Function for collinear equilibrium points
        fun = @(eta)dl_maneuver(eta, delta_i, l);

        % Use fsolve to find root of function
        eta_star(i, j) = fsolve(fun, x0, options);
    end
end

% figure()
% contour(delta_i_range, l_range, eta_star')

%% Plot

figure()
plot(l_range, eta_star(1,:))
hold on
plot(l_range, eta_star(2,:))
plot(l_range, eta_star(3,:))
plot(l_range, eta_star(4,:))
plot(l_range, eta_star(5,:))
plot(l_range, eta_star(6,:),'o')
plot(l_range, eta_star(7,:),'o')
plot(l_range, eta_star(8,:),'o')
plot(l_range, eta_star(9,:),'o')
plot(l_range, eta_star(10,:),'o')
hold off
legend("20°", "40°", "60°", "80°", "100°", "120°", "140°", "160°", "180°")
title("Optimal value of \eta^* vs l with varying \Deltai")
xlabel("l", 'FontSize',15)
ylabel("\eta^*", 'FontSize',15)

%% Function

function out = dl_maneuver(eta, delta_i, l)
    dv1 = sqrt((2*l)/(1+l) + 1 - 2*cos(eta*delta_i)*sqrt((2*l)/(1+l)));
    dv2 = sqrt((4)/(l*(1+l)) - 2*cos((1-eta)*delta_i)/(l*(1+l)));
    out = (delta_i * sqrt((2*l)/(1+l)) * sin(eta*delta_i))/(dv1) - (delta_i*sin((1-eta)*delta_i)*sqrt((2)/(l*(1+l))))/(dv2);
end
    

