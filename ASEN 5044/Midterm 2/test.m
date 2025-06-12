clear; clc; close all;

%% Part c
% Load data
data_part_c = load("midterm2_problem3c.mat");
y_aug_hist = data_part_c.yaugHist;
[p, T] = size(y_aug_hist);

% Set n = 8 because x = [x_a, x_b]^T
n = 8;

% Given
delta_t = 0.5;
omega_a = 0.045;
odt_a = delta_t*omega_a;
omega_b = -0.045;
odt_b = delta_t*omega_b;
R_a_static = [75 7.5; 7.5 75];
R_D = [8000 500; 500 8000];

% Construct F_a, F_b matrices
F_a = [1 sin(odt_a)/omega_a 0 -(1-cos(odt_a))/omega_a;
        0 cos(odt_a) 0 -sin(odt_a);
        0 (1-cos(odt_a))/omega_a 1 sin(odt_a)/omega_a;
        0 sin(odt_a) 0 cos(odt_a)];

F_b = [1 sin(odt_b)/omega_b 0 -(1-cos(odt_b))/omega_b;
        0 cos(odt_b) 0 -sin(odt_b);
        0 (1-cos(odt_b))/omega_b 1 sin(odt_b)/omega_b;
        0 sin(odt_b) 0 cos(odt_b)];

% Big F matrix is [F_a 0; 0 F_b]
F = [F_a, zeros(4,4); zeros(4,4), F_b];

% Construct H matrix
H = [1 zeros(1,7); 0 0 1 zeros(1,5); 1 zeros(1,3) -1 zeros(1,3); 0 0 1 zeros(1,3) -1 0];

for i = 1:T
    % R_a is dynamically changing
    R_a_k((2*i-1):(2*i),:) = R_a_static + [12.5*sin(i/10), 25.5*sin(i/10); 25.5*sin(i/10), 12.5*cos(i/10)];
end

for i = 1:T
    % Full R matrix is [R_a_k, 0; 0, R_D]; 
    R((4*(i-1)+1):(4*i),:) = [R_a_k((2*i-1):(2*i),:), zeros(2,2); zeros(2,2), R_D];
end

% Reshape given y to a column vector
y_vec = reshape(y_aug_hist, [T*p, 1]);

% Initial estimates
x_hat_0_est = zeros(8,1);
P_0_est = eye(n) * 100000;

% Initialize x_hat_0 and P_k matrix
x_hat_0 = [x_hat_0_est];
P_k = [P_0_est];

% Ancilliary identity matrix needed later
I_n = eye(n);

% RLLS loop
for i = 1:T
    % Set k variables
    P_k_minus_1 = P_k((8*(i-1)+1:8*i), :);    
    H_k = H*F^i;
    R_k = R((4*(i-1)+1):(4*i),:);
    x_hat_0_k_minus_1 = x_hat_0(:,i);
    y_k = y_vec(4*(i-1)+1:4*i, :);
    
    % Loop
    K_k = P_k_minus_1 * H_k' * inv(H_k*P_k_minus_1*H_k' + R_k);
    x_hat_0(:, i+1) = x_hat_0_k_minus_1 + K_k*(y_k - H_k*x_hat_0_k_minus_1);
    P_k(8*i+1:8*(i+1),:) = (I_n - K_k*H_k) * P_k_minus_1 * (I_n - K_k*H_k)' + K_k*R_k*K_k';
end

% Get x_hat_0 for aircraft a and b
x_hat_0_a = x_hat_0(1:4,:);
x_hat_0_b = x_hat_0(5:8,:);

% Last columns are the RLLS estimate
x_hat_0_a_rlls = x_hat_0_a(:,end);
x_hat_0_b_rlls = x_hat_0_b(:,end);

%% Part d
% Set k timestep
k = 0:T;

% Plot state estimate for aircraft a
figure()
subplot(2,2,1)
plot(k, x_hat_0_a(1,:))
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,2,2)
plot(k, x_hat_0_a(2,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\xi}$ [m/s]', 'Interpreter','latex')

subplot(2,2,3)
plot(k, x_hat_0_a(3,:))
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")

subplot(2,2,4)
plot(k, x_hat_0_a(4,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\eta}$ [m/s]', 'Interpreter','latex')
sgtitle("Aircraft a state estimate vs k")

% Plot state estimate for aircraft b
figure()
subplot(2,2,1)
plot(k, x_hat_0_b(1,:))
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,2,2)
plot(k, x_hat_0_b(2,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\xi}$ [m/s]', 'Interpreter','latex')

subplot(2,2,3)
plot(k, x_hat_0_b(3,:))
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")

subplot(2,2,4)
plot(k, x_hat_0_b(4,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\eta}$ [m/s]', 'Interpreter','latex')
sgtitle("Aircraft b state estimate vs k")

% Parse through P_k to get 2 sigma values for a and b aircrafts
for i = 1:T+1
    P_k_current = P_k(8*(i-1)+1:8*i,:);
    for j = 1:n
        P_2sig(j) = 2 * sqrt(P_k_current(j,j));
    end
    
    P_2sig_a(:,i) = [P_2sig(1); P_2sig(2); P_2sig(3); P_2sig(4)];
    P_2sig_b(:,i) = [P_2sig(5); P_2sig(6); P_2sig(7); P_2sig(8)];

end

% Plot +2sig bounds for aircraft a
figure()
subplot(2,2,1)
plot(k, P_2sig_a(1,:))
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,2,2)
plot(k, P_2sig_a(2,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\xi}$ [m/s]', 'Interpreter','latex')

subplot(2,2,3)
plot(k, P_2sig_a(3,:))
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")

subplot(2,2,4)
plot(k, P_2sig_a(4,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\eta}$ [m/s]', 'Interpreter','latex')
sgtitle("Aircraft a 2\sigma bounds vs k")

% Plot +2sig bounds for aircraft b
figure()
subplot(2,2,1)
plot(k, P_2sig_b(1,:))
xlabel("k time step [0.5 sec]")
ylabel("\xi [m]")

subplot(2,2,2)
plot(k, P_2sig_b(2,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\xi}$ [m/s]', 'Interpreter','latex')

subplot(2,2,3)
plot(k, P_2sig_b(3,:))
xlabel("k time step [0.5 sec]")
ylabel("\eta [m]")

subplot(2,2,4)
plot(k, P_2sig_b(4,:))
xlabel("k time step [0.5 sec]")
ylabel('$\dot{\eta}$ [m/s]', 'Interpreter','latex')
sgtitle("Aircraft b 2\sigma bounds vs k")
