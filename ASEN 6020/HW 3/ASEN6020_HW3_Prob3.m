clear; clc; close all;

% ASEN 6020 - HW 3 - Problem 3
% Spring 2025
% Jash Bhalavat

A = [0 1; 0 0];
B = [0;1];
S = eye(2);
Q = zeros(2);
R = 1;
x0 = [1;1];

phi_xx = A;
phi_xl = -B*inv(R)*B';
phi_lx = -Q;
phi_ll = -A';

% xf = inv(phi_xx + phi_xl*S)*x0;

KX0 = [reshape(S, [4,1]); x0];

[T, KX] = ode45(@mRiccati, [0 10], KX0, [], A, B, Q);

for i = 1:length(T)
    u_t(i) = -inv(R)*B'*reshape(KX(i,1:4), [2,2])*KX(i,5:6)';
end

function dKXdt = mRiccati(t, KX, A, B, Q)
    K = reshape(KX(1:4), size(A)); % Convert from "n^2"-by-1 to "n"-by-"n"
    X = KX(5:6);
    dKdt =  K*A + A'*K - K*B*B'*K - Q; % Determine derivative
    dXdt = (A - B*B'*K)*X;
    dKdt = dKdt(:); % Convert from "n"-by-"n" to "n^2"-by-1
    dXdt = dXdt(:);
    dKXdt = [dKdt;dXdt];
end

figure()
subplot(2,1,1)
plot(T, u_t)
xlabel("time")
ylabel("u(t)")
grid on

subplot(2,1,2)
plot(T, KX(:,5))
hold on
plot(T, KX(:,6))
hold off
xlabel("time")
ylabel("x(t)")
legend("x1", "x2")
grid on

sgtitle("Optimal Feedback Control Law & Optimal State")
