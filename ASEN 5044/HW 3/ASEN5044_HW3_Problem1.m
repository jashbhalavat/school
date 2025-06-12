clear; clc; close all;

A = [0 1 0 0; -2 0 1 0; 0 0 0 1; 1 0 -2 0];
B = [0 0; -1 0; 0 0; 1 1];

C = [1 0 0 0; 0 1 0 -1];
D = [0 0; 0 0];

H = C;
M = D;

delta_t = 0.05;

% Part a
% A is always square matrix
n = length(A);
m = size(B,2);

A_hat = [A B; zeros([n - m, n+m])];

syms t

stm = expm(A_hat*t);
stmf(t) = stm;
f(t) = stm(1:4, 1:4);
g(t) = stm(1:4, 5:6);


% Part b
% is u = 0?

F = double(f(delta_t));
G = double(g(delta_t));

O = [H; H*F; H*F*F; H*F*F*F];

rank_O = rank(O);

% Part c
data = load("hw3problem1data.mat");
U = data.Udata;
Y = data.Ydata;

t0 = 0;
tf = 5;

time_steps = t0:delta_t:tf;

% Estimating x(0)
Y_temp = Y';
Y_sole = reshape(Y_temp, [2*length(Y_temp), 1]);

U_temp = U(1:end-1,:)';
U_sole = reshape(U_temp, [2*length(U_temp), 1]);

n = 3;

O_sole = [];
for i = 1:n
    O_sole = [O_sole; H*F^i];
end

p = 2;
HFG_sole = [];
for i = 1:n
    temp_vec = zeros([p, 2*n]);
    count = 1;
    for j = 1:i
        temp_vec(:,count:count+1) = H*F^(i-j)*G;
        count = count + 2;
    end
    HFG_sole = [HFG_sole; temp_vec];
end

LHS = Y_sole(1:2*n) - HFG_sole*U_sole(1:2*n);
gram = O_sole' * O_sole;
x0 = inv(gram) * O_sole' * LHS;

x = x0;

for i = 1:length(time_steps)
    x(:,i+1) = F * x(:,i) + G * U(i,:)';
    y_predicted(:,i) = H * x(:,i) + M * U(i,:)';
end

figure()
plot(time_steps(2:end), x(1,2:length(time_steps)))
hold on
plot(time_steps(2:end), x(2,2:length(time_steps)))
plot(time_steps(2:end), x(3,2:length(time_steps)))
plot(time_steps(2:end), x(4,2:length(time_steps)))
hold off
legend("q1 [m]", "q1dot [m/s]", "q2 [m]", "q2dot[m/s]")
xlabel("Time [seconds]")
title("Estimated x states")

figure()
plot(time_steps(2:end), y_predicted(1,2:end))
hold on
plot(time_steps(2:end), y_predicted(2,2:end))
hold off
legend("y_1 [m]", "y_2 [m/s]")
xlabel("Time [seconds]")
title("Predicted y states")

figure()
subplot(2,1,1)
plot(time_steps(2:end), Y(:,1)-y_predicted(1,2:end)')
xlabel("Time [seconds]")
ylabel("Y1 [m]")
title("Predict y1 - recorded y1")

subplot(2,1,2)
plot(time_steps(2:end), Y(:,2)-y_predicted(2,2:end)')
xlabel("Time [seconds]")
ylabel("Y2 [m/s]")
title("Predict y2 - recorded y2")


% Part e
H_single = [1 0 0 0; 1 0 0 0; 1 0 0 0];
O_single = [H_single*F; H_single*F*F; H_single*F*F*F; H_single*F*F*F*F];

rank_O_single = rank(O_single);





