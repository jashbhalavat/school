clear; clc; close all;

% ASEN 5044, HW 7
% Fall 2024
% Jash Bhalavat

syms delta_t
syms omega

A = [0 1 0 0; 0 0 0 -omega; 0 0 0 1; 0 omega 0 0];

e_at = expm(A * delta_t);

delta_t = 0.5;
omega = 0.045;

% Omega times delta_t
odt = omega * delta_t;

e_at_given = [1 sin(odt)/omega 0 -(1 - cos(odt))/omega; 0 cos(odt) 0 -sin(odt); 0 (1-cos(odt))/omega 1 sin(odt)/omega; 0 sin(odt) 0 cos(odt)];

u_0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P_a_0 = diag([10, 2, 10, 2]);

T = 300;

F = e_at_given;

mu = [u_0];
P = [P_a_0];
P_1_2sig = 2 * sqrt(P(1, 1));
P_2_2sig = 2 * sqrt(P(2, 2));
P_3_2sig = 2 * sqrt(P(3, 3));
P_4_2sig = 2 * sqrt(P(4, 4));

for i = 1:T
    mu(:, i+1) = F^i * u_0;
    F_t = F';
    P_temp = F^i * P_a_0 * F_t^i;
    P(:, i+1) = get_diag(P_temp);
    P_1_2sig(i+1) = 2 * sqrt(P(1, i+1));
    P_2_2sig(i+1) = 2 * sqrt(P(2, i+1));
    P_3_2sig(i+1) = 2 * sqrt(P(3, i+1));
    P_4_2sig(i+1) = 2 * sqrt(P(4, i+1));
end

function out = get_diag(mat)
    len = length(mat);
    for i = 1:len
        out(i, 1) = mat(i, i);
    end
end


k = 0:T;

figure()
plot(k, mu(1, :), 'r')
hold on
plot(k, mu(1, :) + P_1_2sig, '--k')
plot(k, mu(1, :) - P_1_2sig, '--k')
hold off
title("Predicted Results for \DeltaT = 0.5 sec")
legend("mean \xi(k)", "2\sigma bounds", 'FontSize', 11)
xlabel("DT Time Steps K [0.5 sec]")
ylabel("\xi_k [m]", 'FontSize', 15)


figure()
plot(k, mu(2, :), 'r')
hold on
plot(k, mu(2, :) + P_2_2sig, '--k')
plot(k, mu(2, :) - P_2_2sig, '--k')
hold off
title("Predicted Results for \DeltaT = 0.5 sec")
legend("mean dot \xi(k)", "2\sigma bounds", 'FontSize', 11)
xlabel("DT Time Steps K [0.5 sec]")
ylabel('$\dot{\xi_k}$ [m/s]', 'Interpreter','latex', 'FontSize', 15)


figure()
plot(k, mu(3, :), 'r')
hold on
plot(k, mu(3, :) + P_3_2sig, '--k')
plot(k, mu(3, :) - P_3_2sig, '--k')
hold off
title("Predicted Results for \DeltaT = 0.5 sec")
legend("mean \eta(k)", "2\sigma bounds", 'FontSize', 11)
xlabel("DT Time Steps K [0.5 sec]")
ylabel("\eta_k [m]", 'FontSize', 15)


figure()
plot(k, mu(4, :), 'r')
hold on
plot(k, mu(4, :) + P_4_2sig, '--k')
plot(k, mu(4, :) - P_4_2sig, '--k')
hold off
title("Predicted Results for \DeltaT = 0.5 sec")
legend("mean dot \eta(k)", "2\sigma bounds", 'FontSize', 11)
xlabel("DT Time Steps K [0.5 sec]")
ylabel('$\dot{\eta_k}$ [m/s]', 'Interpreter','latex', 'FontSize', 15)


figure()
plot(k, P_1_2sig)
title("Positive 2\sigma vs Time")
xlabel("DT Time Steps K [0.5 sec]")
ylabel("\xi_k [m]", 'FontSize', 15)

figure()
plot(k, P_2_2sig)
title("Positive 2\sigma vs Time")
xlabel("DT Time Steps K [0.5 sec]")
ylabel('$\dot{\xi_k}$ [m/s]', 'Interpreter','latex', 'FontSize', 15)
ylim([0 2*P_2_2sig(1)])

figure()
plot(k, P_3_2sig)
title("Positive 2\sigma vs Time")
xlabel("DT Time Steps K [0.5 sec]")
ylabel("\eta_k [m]", 'FontSize', 15)

figure()
plot(k, P_4_2sig)
title("Positive 2\sigma vs Time")
xlabel("DT Time Steps K [0.5 sec]")
ylabel('$\dot{\eta_k}$ [m/s]', 'Interpreter','latex', 'FontSize', 15)
ylim([0 2*P_4_2sig(1)])

%% Problem 3

delta_t = 0.5;

omega_a = 0.045;
omega_b = -0.045;

odt_a = omega_a * delta_t;

F_a = [1 sin(odt_a)/omega_a 0 -(1-cos(odt_a))/omega_a;
        0 cos(odt_a) 0 -sin(odt_a);
        0 (1-cos(odt_a))/omega_a 1 sin(odt_a)/omega_a;
        0 sin(odt_a) 0 cos(odt_a)];

odt_b = omega_b * delta_t;

F_b = [1 sin(odt_b)/omega_b 0 -(1-cos(odt_b))/omega_b;
        0 cos(odt_b) 0 -sin(odt_b);
        0 (1-cos(odt_b))/omega_b 1 sin(odt_b)/omega_b;
        0 sin(odt_b) 0 cos(odt_b)];

u_a_0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
u_b_0 = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)];

u_a = [u_a_0];
u_b = [u_b_0];

P_a_0 = diag([10, 4, 10, 4]);
P_b_0 = diag([11, 3.5, 11, 3.5]);

P_a = [P_a_0];
P_b = [P_b_0];

T = 150;
t = 0:delta_t:T;

for i = 1:T/delta_t
    u_a(:, i+1) = F_a^i * u_a_0;
    P_a((4*i + 1):(4*(i+1)), :) = F_a^i * P_a_0 * (F_a')^i;

    u_b(:, i+1) = F_b^i * u_b_0;
    P_b((4*i + 1):(4*(i+1)), :) = F_b^i * P_b_0 * (F_b')^i;
end

u_c = u_a - u_b;
P_c = P_a + P_b;

u_r_c = [u_c(1, :); u_c(3, :)];
P_r_c = [P_c(1, 1), P_c(1, 3); P_c(3, 1), P_c(3, 3)];

xi_R = 100;
eta_R = 100;
x_R = [xi_R; eta_R];

cdf_x_R = mvncdf(x_R, u_r_c(:,1), P_r_c);
cdf_neg_x_R = mvncdf(-x_R, u_r_c(:,1), P_r_c);

probability_of_collision = cdf_x_R - cdf_neg_x_R;

for i = 1:T/delta_t
    P_c_temp = [P_c((4*i + 1):(4*(i+1)),:)];
    P_r_c_current = [P_c_temp(1,1), P_c_temp(1,3); P_c_temp(3,1), P_c_temp(3,3)];
    P_r_c((2*i + 1):(2*(i+1)), :) = P_r_c_current;

    cdf_x_R = mvncdf(x_R, u_r_c(:,i), P_r_c_current);
    cdf_neg_x_R = mvncdf(-x_R, u_r_c(:,i), P_r_c_current);

    probability_of_collision(i+1) = cdf_x_R - cdf_neg_x_R;
end

figure()
plot(t, probability_of_collision*100, 'LineWidth',1.5)
xlabel("Time [sec]")
ylabel("Probability [%]")
ylim([0, 100])
title("Probability of collision of the two aircrafts")







