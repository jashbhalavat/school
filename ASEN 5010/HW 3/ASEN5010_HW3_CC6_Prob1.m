clear; clc; close all;

% ASEN 5010 - CC6 Parallel Axis Theorem, Problem 1
% Spring 2025
% Jash Bhalavat

D2R = pi/180;
t1 = -10*D2R;
t2 = 10*D2R;
t3 = 5*D2R;

BN = R1(t3)*R2(t2)*R3(t1);

I_C_B = [10, 1, -1;
    1, 5, 1;
    -1, 1, 8];

mass = 12.5;

R_CP_N = [-0.5, 0.5, 0.25];
R_CP_B = BN * R_CP_N';

I_P_B = I_C_B + mass * skew_symmetric(R_CP_B) * skew_symmetric(R_CP_B)';

% Ip_B = [12.3212520668661, 4.19755561866020, -0.158138667715023,
              % 4.19755561866020, 9.86047156663829, 0.428471419436600,
              % -0.158138667715023, 0.428471419436600, 14.8807763664957] # kg m^2


%%
clear; clc; close all;

% ASEN 5010 - CC6.1 Coordinate Transformation
% Spring 2025
% Jash Bhalavat

%% Problem 1

I_C_B = [10, 1, -1;
    1, 5, 1;
    -1, 1, 8];
sigma_DB = [0.1, 0.2, 0.3];
sigma_DB_sq = norm(sigma_DB)^2;
beta_DB_0 = (1-sigma_DB_sq)/(1+sigma_DB_sq);
beta_DB_i = 2/(1+sigma_DB_sq) * sigma_DB;
beta_DB = [beta_DB_0; beta_DB_i'];
DB = ep2dcm(beta_DB);

I_C_D = DB * I_C_B * DB';

% Ic_D = [5.42779505231195, -1.77341011998767, 1.37988230580880,
              % -1.77341011998767, 9.27952214100775, -0.530473519280645,
              % 1.37988230580880, -0.530473519280645, 8.29268280668029] # kg m^2


function out = ep2dcm(beta)
    beta0 = beta(1);
    beta1 = beta(2);
    beta2 = beta(3);
    beta3 = beta(4);
    out = [beta0^2 + beta1^2 - beta2^2 - beta3^2, 2*(beta1*beta2 + beta0*beta3), 2*(beta1*beta3 - beta0*beta2);
        2*(beta1*beta2 - beta0*beta3), beta0^2 - beta1^2 + beta2^2 - beta3^2, 2*(beta2*beta3 + beta0*beta1);
        2*(beta1*beta3 + beta0*beta2), 2*(beta2*beta3 - beta0*beta1), beta0^2 - beta1^2 - beta2^2 + beta3^2];
end

%% Problem 2
eigvals = eig(I_C_B);
% principalInertias = [10.4741936586104, 8.11268085340805, 4.41312548798157 ] # kg m^2

%% Problem 3
[V, D] = eig(I_C_B);
FB = [V(:,3)';
        V(:,2)';
        V(:,1)'];
test = FB * I_C_B * FB';

% FB = [0.936164162496731, 0.110017824614094, -0.333905284660151,
              % 0.272608605545180, 0.372563629308523, 0.887063070079673,
              % 0.221993713963946, -0.921462110118286, 0.318788912255193]
