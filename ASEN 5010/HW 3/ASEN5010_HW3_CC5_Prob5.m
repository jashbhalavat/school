clear; clc; close all;

% ASEN 5010 - HW 3, Concept Check 5 Quest, Problem 5
% Spring 2025
% Jash Bhalavat

v1_B = [0.8273; 0.5541; -0.0920];
v1_B_hat = normalize(v1_B);
v2_B = [-0.8285; 0.5522; -0.0955];
v2_B_hat = normalize(v2_B);

v1_N = [-0.1517; -0.9669; 0.2050];
v1_N_hat = normalize(v1_N);
v2_N = [-0.8393; 0.4494; -0.3044];
v2_N_hat = normalize(v2_N);

B = v1_B_hat * v1_N_hat' + v2_B_hat * v2_N_hat';
Z = [B(2,3) - B(3,2); B(3,1) - B(1,3); B(1,2) - B(2,1)];
S = B + B';
sigma = trace(B);

K = [sigma Z';
    Z, S - sigma*eye(3)];

lambda_opt = 2;

q_bar = inv((lambda_opt + sigma)*eye(3)) * Z;
beta_bar = (1/(sqrt(1+q_bar'*q_bar))) * [1; q_bar];

BN_est = ep2dcm(beta_bar);

%% NR method

syms f(s)

f(s) = det(K - s*eye(4));
df = diff(f, s);

tol = 1e-12;
counter = 1;
lambda = lambda_opt;
error = 1;

while ((abs(error) > tol) && (counter < 10))
    lambda(counter + 1) = lambda(counter) - double(f(lambda(counter)))/double(df(lambda(counter)));
    error = f(lambda(counter+1));
    counter = counter + 1;
end

q_bar = inv((lambda(end) + sigma)*eye(3) - S) * Z;
beta_bar = (1/(sqrt(1+q_bar'*q_bar))) * [1; q_bar];

BN_est = ep2dcm(beta_bar);

% dcm = [0.827300105095276,-0.485552816844880, 0.282511943399475,
%        0.508424652181109, 0.861059924503972, -0.00895429884087926,
       % -0.238911927629962, 0.151043928932001, 0.959221988055382]

% 2nd attempt
% dcm = [0.415936375258194, -0.854893534518689, 0.310087046448601,
              % -0.833756672976140, -0.494636774164877, -0.245324829382438,
              % 0.363107066861598, -0.156497623885817, -0.918510616004660]

function hat = normalize(vec)
    vec_norm = norm(vec);
    hat = vec ./ vec_norm;
end

function out = ep2dcm(beta)
    beta0 = beta(1);
    beta1 = beta(2);
    beta2 = beta(3);
    beta3 = beta(4);
    out = [beta0^2 + beta1^2 - beta2^2 - beta3^2, 2*(beta1*beta2 + beta0*beta3), 2*(beta1*beta3 - beta0*beta2);
        2*(beta1*beta2 - beta0*beta3), beta0^2 - beta1^2 + beta2^2 - beta3^2, 2*(beta2*beta3 + beta0*beta1);
        2*(beta1*beta3 + beta0*beta2), 2*(beta2*beta3 - beta0*beta1), beta0^2 - beta1^2 - beta2^2 + beta3^2];
end