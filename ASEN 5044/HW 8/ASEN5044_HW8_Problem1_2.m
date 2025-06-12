clear; clc; close all;

% ASEN 5044 - HW 8
% Fall 2024
% Jash Bhalavat

%% Problem 1
% Given
delta_t = 0.5;
omega_a = 0.045;
odt_a = delta_t*omega_a;
omega_b = -0.045;
odt_b = delta_t*omega_b;

A_a = [0 1 0 0; 0 0 0 -omega_a; 0 0 0 1; 0 omega_a 0 0];
A_b = [0 1 0 0; 0 0 0 -omega_b; 0 0 0 1; 0 omega_b 0 0];

n = length(A_a);

% Construct F_a, F_b matrices
F_a = [1 sin(odt_a)/omega_a 0 -(1-cos(odt_a))/omega_a;
        0 cos(odt_a) 0 -sin(odt_a);
        0 (1-cos(odt_a))/omega_a 1 sin(odt_a)/omega_a;
        0 sin(odt_a) 0 cos(odt_a)];

F_b = [1 sin(odt_b)/omega_b 0 -(1-cos(odt_b))/omega_b;
        0 cos(odt_b) 0 -sin(odt_b);
        0 (1-cos(odt_b))/omega_b 1 sin(odt_b)/omega_b;
        0 sin(odt_b) 0 cos(odt_b)];

q_omega = 10;
W = q_omega*[2 0.05; 0.05 0.5];

gamma_a = [0 0; 1 0; 0 0; 0 1];
gamma_b = [0 0; 1 0; 0 0; 0 1];

Z_a = delta_t * [-A_a gamma_a*W*gamma_a'; zeros(n), A_a'];
Z_b = delta_t * [-A_b gamma_b*W*gamma_b'; zeros(n), A_b'];

e_z_a = expm(Z_a);
e_z_b = expm(Z_b);

F_inv_Q_a = e_z_a(1:4, 5:8);
F_inv_Q_b = e_z_b(1:4, 5:8);
F_a_t = e_z_a(5:8, 5:8);
F_b_t = e_z_b(5:8, 5:8);

Q_a = F_a_t' * F_inv_Q_a;
Q_b = F_b_t' * F_inv_Q_b;


%% Problem 2

rng(100);
H = [1 0 0 0; 0 0 1 0];
R_a = [20 0.05; 0.05 20];

data = load("hw8problemdata.mat");
x_a_single_truth = data.xasingle_truth;

p = size(H,1);
% Subtracting 1 because x_a_single_truth starts at 0
T = size(x_a_single_truth, 2) - 1;

% Part a
S_v_a = chol(R_a, 'lower');

% Necessary variables
I_p = eye(p);
zeros_p = zeros(p,1);

for i = 1:T
    q_k_a = mvnrnd(zeros_p, I_p)';
    % Using x(:,i+1) because x starts at 0
    y_a_k(:,i) = H*x_a_single_truth(:,i+1) + S_v_a*q_k_a;
end

k_20_sec = 1:40;

figure()
subplot(2,1,1)
plot(k_20_sec, y_a_k(1,1:40))
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")
% hold on
% plot(time_20_sec, x_a_single_truth(1,2:41))
% hold off

subplot(2,1,2)
plot(k_20_sec, y_a_k(2,1:40))
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")
% hold on
% plot(time_20_sec, x_a_single_truth(3,2:41))
% hold off
sgtitle("Simulated noisy measurements y_a(k)")

%% Part b
mu_a_0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P_a_0 = 900 * diag([10, 2, 10, 2]);

x_a_k = mu_a_0;
P_a_k = P_a_0;

Qkf_a = Q_a;
Rkf_a = R_a;

tvec = 1:T;

G_a = zeros(4,1);
u_a = zeros(1,T);

[x_a_kf, P_a_kf] = kalman_filter_hw8(tvec, F_a, G_a, x_a_k, u_a, P_a_k, Qkf_a, Rkf_a, y_a_k, H);

figure()
subplot(2,2,1)
plot(tvec, x_a_single_truth(1,2:end)-x_a_kf(1,:))
hold on
plot(tvec,2*sqrt(squeeze(P_a_kf(1,1,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_a_kf(1,1,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,2,2)
plot(tvec, x_a_single_truth(2,2:end)-x_a_kf(2,:))
hold on
plot(tvec,2*sqrt(squeeze(P_a_kf(2,2,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_a_kf(2,2,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\xi}$ [m/s]', 'Interpreter','latex')
legend("Error", "2\sigma bounds")

subplot(2,2,3)
plot(tvec, x_a_single_truth(3,2:end)-x_a_kf(3,:))
hold on
plot(tvec,2*sqrt(squeeze(P_a_kf(3,3,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_a_kf(3,3,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")

subplot(2,2,4)
plot(tvec, x_a_single_truth(4,2:end)-x_a_kf(4,:))
hold on
plot(tvec,2*sqrt(squeeze(P_a_kf(4,4,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_a_kf(4,4,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\eta}$ [m/s]', 'Interpreter','latex')
sgtitle("Aircraft A Estimated State Error and Estimated 2\sigma Bounds")




