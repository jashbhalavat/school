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

n = length(Y);

O_sole = [];
for i = 1:length(Y)
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

LHS = Y_sole - HFG_sole*U_sole;
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
legend("q1", "q1dot", "q2", "q2dot")
xlabel("Time [seconds]")
title("Estimated x states")

figure()
plot(time_steps(2:end), Y(:,1))
hold on
plot(time_steps(2:end), Y(:,2))
plot(time_steps(2:end), y_predicted(1,2:end))
plot(time_steps(2:end), y_predicted(2,2:end))
hold off
legend("y1", "y2", "y_predicted1", "y_predicted2")
xlabel("Time [seconds]")
title("Predicted y states")


