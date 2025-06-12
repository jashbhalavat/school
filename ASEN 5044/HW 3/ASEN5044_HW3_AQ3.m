clear; clc; close all;

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

F = double(f(delta_t));
G = double(g(delta_t));

cap_C = [G, F*G, F*F*G, F*F*F*G];
rank_cap_c = rank(cap_C);

% Part b
G_u2_only = G(:,2)
cap_C_u2_only = [G_u2_only, F*G_u2_only, F*F*G_u2_only, F*F*F*G_u2_only];
rank_cap_c_u2_only = rank(cap_C_u2_only);

G_u1_only = G(:,1)
cap_C_u1_only = [G_u1_only, F*G_u1_only, F*F*G_u1_only, F*F*F*G_u1_only];
rank_cap_c_u1_only = rank(cap_C_u1_only);


