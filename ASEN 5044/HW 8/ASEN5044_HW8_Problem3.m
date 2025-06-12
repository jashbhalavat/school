clear; clc; close all;
rng(100);

% ASEN 5044 - HW 8 Problem 3
% Fall 2024, Jash Bhalavat

% From problem 1
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

%% Problem 3 Part a

data = load("hw8problemdata.mat");
x_a_double_truth = data.xadouble_truth;
x_b_double_truth = data.xbdouble_truth;

% Simulate noisy measurements for a'
H_a = [1 0 0 0; 0 0 1 0];
R_a = [20 0.05; 0.05 20];

p = size(H_a,1);
% Subtracting 1 because x_a_single_truth starts at 0
T = size(x_a_double_truth, 2) - 1;

% Part a
S_v_a = chol(R_a, 'lower');

% Necessary variables
I_p = eye(p);
zeros_p = zeros(p,1);

for i = 1:T
    q_k_a = mvnrnd(zeros_p, I_p)';
    % Using x(:,i+1) because x starts at 0
    y_a_k(:,i) = H_a*x_a_double_truth(:,i+1) + S_v_a*q_k_a;
end

% Simulate y_d noisy measurements
x_truth = [x_a_double_truth; x_b_double_truth];
H_d = [1 0 0 0 -1 0 0 0; 0 0 1 0 0 0 -1 0];
R_d = [10 0.15; 0.15 10];

S_v_d = chol(R_d, 'lower');

for i = 1:T
    q_k_d = mvnrnd(zeros_p, I_p)';
    y_d_k(:,i) = H_d*x_truth(:,i+1) + S_v_d*q_k_d;
end

y_s = [y_a_k; y_d_k];

tvec = 1:T;

% figure()
% plot(tvec, y_d_k(1,:))
% hold on
% plot(tvec, x_a_double_truth(1,2:end)-x_b_double_truth(1,2:end))
% hold off

mu_a_0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P_a_0 = 900 * diag([10, 2, 10, 2]);
mu_b_0 = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)];
P_b_0 = 900 * diag([11, 4, 11, 4]);

F = [F_a, zeros(4,4); zeros(4,4), F_b];
G = zeros(8,1);
u = zeros(1,T);

xk = [mu_a_0; mu_b_0];
Pk = [P_a_0 zeros(4,4); zeros(4,4) P_b_0];

Qkf = [Q_a, zeros(4,4); zeros(4,4) Q_b];
Rkf = [R_a, zeros(2,2); zeros(2,2) R_d];

H_s = [H_a, zeros(2,4); H_d];

[x_kf, P_kf] = kalman_filter_hw8(tvec, F, G, xk, u, Pk, Qkf, Rkf, y_s, H_s);

figure()
subplot(2,1,1)
plot(tvec, x_a_double_truth(1,2:end)-x_kf(1,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf(1,1,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf(1,1,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,1,2)
plot(tvec, x_a_double_truth(3,2:end)-x_kf(3,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf(3,3,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf(3,3,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")
sgtitle("Aircraft A position errors and 2\sigma bounds")
legend("Error", "2\sigma bounds")

figure()
subplot(2,1,1)
plot(tvec, x_b_double_truth(1,2:end)-x_kf(5,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf(5,5,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf(5,5,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,1,2)
plot(tvec, x_b_double_truth(3,2:end)-x_kf(7,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf(7,7,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf(7,7,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")
sgtitle("Aircraft B position errors and 2\sigma bounds")
legend("Error", "2\sigma bounds")

%% Problem 3 Part b

[x_kf_partb, P_kf_partb] = kalman_filter_hw8(tvec, F, G, xk, u, Pk, Qkf, R_d, y_d_k, H_d);

figure()
subplot(2,1,1)
plot(tvec, x_a_double_truth(1,2:end)-x_kf_partb(1,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf_partb(1,1,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf_partb(1,1,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,1,2)
plot(tvec, x_a_double_truth(3,2:end)-x_kf_partb(3,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf_partb(3,3,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf_partb(3,3,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")
sgtitle("Aircraft A position errors and 2\sigma bounds")
legend("Error", "2\sigma bounds")

figure()
subplot(2,1,1)
plot(tvec, x_b_double_truth(1,2:end)-x_kf_partb(5,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf_partb(5,5,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf_partb(5,5,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,1,2)
plot(tvec, x_b_double_truth(3,2:end)-x_kf_partb(7,:))
hold on
plot(tvec,2*sqrt(squeeze(P_kf_partb(7,7,:))'),'b--','LineWidth',1.25)
plot(tvec,-2*sqrt(squeeze(P_kf_partb(7,7,:))'),'b--','LineWidth',1.25)
hold off
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")
sgtitle("Aircraft B position errors and 2\sigma bounds")
legend("Error", "2\sigma bounds")

