BN = [1 0 0; 0 cos(pi/2) sin(pi/2); 0 -sin(pi/2) cos(pi/2)]
RN = [cos(-pi/2) 0 -sin(-pi/2); 0 1 0; sin(-pi/2) 0 cos(-pi/2)]
BN * RN'

clear; clc;

b1 = 1/3 * [ 1 2 -2]
b2 = 1/sqrt(2) * [0, 1, 1]
b3 = 1/(3*sqrt(2)) * [4 -1 1]
BN = [b1; b2; b3]
f1 = 1/4 * [3, -2, sqrt(3)]
f2 = 1/2 * [-1 0 sqrt(3)]
f3 = -1/4 * [sqrt(3), 2*sqrt(3), 1]
FN = [f1; f2; f3]
BN * FN'

dcm = [0, 0, -1,
           1, 0, 0,
           0, -1, 0]
BF = [-0.3720, -0.7440, -0.5550,
              -0.0474, 0.6124, -0.7891,
              0.9270, -0.2673, -0.2630]

BN = [0.3333, 0.6667, -0.6667,
           0, 0.7071, 0.7071,
           0.9428, -0.2357, 0.2357]



clear; clc

BN = [-0.87097, 0.45161, 0.19355; -0.19355, -0.67742, 0.70968; 0.45161, 0.58065, 0.67742];
omega_BN = [0.1; 0.2; 0.3];

omega_BN_tilde = [0 -omega_BN(3) omega_BN(2); omega_BN(3), 0, -omega_BN(1); -omega_BN(2) omega_BN(1), 0];
BN_dot = -omega_BN_tilde * BN;

clear; clc
% CC5, 6
beta = [0.235702, 0.471405, -0.471405, 0.707107];
beta2dcm = ep2dcm(beta);

function out = ep2dcm(beta)
    beta0 = beta(1);
    beta1 = beta(2);
    beta2 = beta(3);
    beta3 = beta(4);
    out = [beta0^2 + beta1^2 - beta2^2 - beta3^2, 2*(beta1*beta2 + beta0*beta3), 2*(beta1*beta3 - beta0*beta2);
        2*(beta1*beta2 - beta0*beta3), beta0^2 - beta1^2 + beta2^2 - beta3^2, 2*(beta2*beta3 + beta0*beta1);
        2*(beta1*beta3 + beta0*beta2), 2*(beta2*beta3 - beta0*beta1), beta0^2 - beta1^2 - beta2^2 + beta3^2];
end

% dcm = [-0.4444, -0.1111, 0.8889,
%               -0.7778, -0.4444, -0.4444,
%               0.4444, -0.8889, 0.1111]

clear; clc
BN = [-0.529403, -0.467056, 0.708231;
    -0.474115, -0.529403, -0.703525;
    0.703525, -0.708231, 0.0588291];
betafromdcm = dcm2ep(BN);
% ep = [0.0024, 0.4896, -0.4896, 0.7344 ];
beta_sm = sheppards_method(BN);
% ep = [0.0024, 0.4851, -0.4851, 0.7276 ];

function out = sheppards_method(C)
    beta0_sq = 1/4 * (1 + trace(C));
    beta1_sq = 1/4 * (1 + 2*C(1,1) - trace(C));
    beta2_sq = 1/4 * (1 + 2*C(2,2) - trace(C));
    beta3_sq = 1/4 * (1 + 2*C(3,3) - trace(C));
    [val, max_idx] = max([beta0_sq, beta1_sq, beta2_sq, beta3_sq]);
    if max_idx == 1
        beta0 = sqrt(beta0_sq);
        beta1 = (C(2,3) - C(3,2))/(4*beta0);
        beta2 = (C(3,1) - C(1,3))/(4*beta0);
        beta3 = (C(1,2) - C(2,1))/(4*beta0);
    elseif max_idx == 2
        beta1 = sqrt(beta1_sq);
        beta0 = (C(2,3) - C(3,2))/(4*beta1);
        beta2 = (C(1,2) + C(2,1))/(4*beta1);
        beta3 = (C(3,1) + C(1,3))/(4*beta1);
    elseif max_idx == 3
        beta2 = sqrt(beta2_sq);
        beta0 = (C(3,1) - C(1,3))/(4*beta2);
        beta1 = (C(1,2) + C(2,1))/(4*beta2);
        beta3 = (C(2,3) + C(3,2))/(4*beta2);
    elseif max_idx == 4
        beta3 = sqrt(beta3_sq);
        beta0 = (C(1,2) - C(2,1))/(4*beta3);
        beta1 = (C(3,1) + C(1,3))/(4*beta3);
        beta2 = (C(2,3) + C(3,2))/(4*beta3);
    end
    out = [beta0, beta1, beta2, beta3];
    if beta0 < 0
        out = -out;
    end
end

function out = dcm2ep(C)
    beta0 = 1/2 * sqrt(C(1,1) + C(2,2) + C(3,3) + 1);
    beta1 = (C(2,3) - C(3,2))/(4*beta0);
    beta2 = (C(3,1) - C(1,3))/(4*beta0);
    beta3 = (C(1,2) - C(2,1))/(4*beta0);
    out = [beta0, beta1, beta2, beta3];
end

clear; clc
D2R = pi/180;
euler = euler_rotation_321(20*D2R, 10*D2R, -10*D2R);
beta = dcm2ep(euler);

% eulerParameters = [0.9760, -0.1006, 0.0704, 0.1798 ]

clear; clc
% CC7 EP
beta_BN = [0.774597,0.258199,0.516398,0.258199];
beta_FB = [0.359211,0.898027,0.179605,0.179605];
dcm_BN = ep2dcm(beta_BN);
dcm_FB = ep2dcm(beta_FB);
dcm_FN = dcm_FB * dcm_BN;
beta_FN = sheppards_method(dcm_FN);
% beta_FN = beta_add_mat(beta_FB) * beta_BN';
% beta_FN = beta_add_mat(beta_BN, beta_FB);
% beta_FN = [-0.0927, 0.8347, 0.5101, -0.1855 ] 
% beta_FN = [0.0927, -0.8347, -0.5101, 0.1855 ] 

function out = beta_add_mat(beta_p, beta_pp)
    mat = [beta_pp(1), -beta_pp(2), -beta_pp(3), -beta_pp(4);
            beta_pp(2), beta_pp(1), beta_pp(4), -beta_pp(3);
            beta_pp(3), -beta_pp(4), beta_pp(1), beta_pp(2)
            beta_pp(4), beta_pp(3), -beta_pp(2), beta_pp(1)];
    out = mat * beta_p';
end

clear; clc;
beta_FN = [0.359211,0.898027,0.179605,0.179605];
beta_BN = [-0.377964,0.755929,0.377964,0.377964];
dcm_FN = ep2dcm(beta_FN);
dcm_BN = ep2dcm(beta_BN);
dcm_FB = dcm_FN * dcm_BN';
beta_FB = sheppards_method(dcm_FB);
% beta_FB = [0.6788, -0.6110, -0.4073, 0 ] 

clear; clc;
% CC8 - EP



clear; clc
% CC7
D2R = pi/180;
dcm = euler_rotation_321(20*D2R, 10*D2R, -10*D2R);
Omega = atan2(dcm(3,1), -dcm(3,2));
i = acos(dcm(3,3));
omega = atan2(dcm(1,3), dcm(2,3));

% eulerAngles = [2.7129, 0.2462, -2.3485 ] # radians

% Test
dcm = euler_rotation_321(60*D2R, 50*D2R, 70*D2R);
Omega = atan2d(dcm(3,1), -dcm(3,2));
i = acosd(dcm(3,3));
omega = atan2d(dcm(1,3), dcm(2,3));


clear; clc
% CC8
D2R = pi/180;
dcm = euler_rotation_321(10*D2R, 20*D2R, 30*D2R);
Omega = atan2d(dcm(3,1), -dcm(3,2));
i = acosd(dcm(3,3));
omega = atan2d(dcm(1,3),dcm(2,3));
% eulerAngles = [40.6423, 35.5313, -36.0524 ] # degrees


dcm_BN = euler_rotation_321(10*D2R, 20*D2R, 30*D2R);
dcm_RN = euler_rotation_321(-5*D2R, 5*D2R, 5*D2R);
dcm_BR = dcm_BN * dcm_RN';
t1 = atan(dcm_BR(1,2)/dcm_BR(1,1))/D2R;
t2 = -asin(dcm_BR(1,3))/D2R;
t3 = atan(dcm_BR(2,3)/dcm_BR(3,3))/D2R;
% eulerAngles = [13.2238, 16.3683, 23.6176 ] # degrees

