clear; clc; close all;

% ASEN 5010 - HW 3, Concept Checks
% Spring 2025
% Jash Bhalavat

%% CC11
crp = [0.1; 0.2; 0.3];
dcm = crp2dcm(crp);
%     dcm = [0.771929824561404,	0.561403508771930,	-0.298245614035088,
       % -0.491228070175439,	0.824561403508772,	0.280701754385965,
       % 0.403508771929825,	-0.0701754385964912,	0.912280701754386]

function dcm = crp2dcm(crp)
    dcm = 1/(1 + crp'*crp) * ((1-crp'*crp)*eye(3) + 2 * crp * crp' -2 * skew_symmetric(crp));
end

function out = skew_symmetric(vec)
    out = [0 -vec(3) vec(2);
            vec(3) 0 -vec(1);
            -vec(2) vec(1) 0];
end

clear; clc;
dcm_bn = [0.333333, 0.871795, -0.358974;
    -0.666667, 0.487179, 0.564103;
    0.666667, 0.0512821, 0.74359]';
ep = sheppards_method(dcm_bn);
crp_bn = [ep(2)/ep(1), ep(3)/ep(1), ep(4)/ep(1)];
% crp = [-0.2, -0.4, -0.6 ]

% 3
% crp = [-0.1, -0.2, -0.3 ]

% 4
clear; clc;
q_fn = [0.1; 0.2; 0.3];
q_bn = [-0.3; 0.3; 0.1];
dcm_fn = crp2dcm(q_fn);
dcm_bn = crp2dcm(q_bn);
dcm_bf = dcm_bn * dcm_fn';
ep_bf = sheppards_method(dcm_bf);
crp_bf = [ep_bf(2)/ep_bf(1), ep_bf(3)/ep_bf(1), ep_bf(4)/ep_bf(1)];
% crp = [-0.311320754716981, 0.188679245283019, -0.273584905660377]
% q_bf = q_bn + -q_fn
q_p = q_bn;
q_pp = -q_fn;
q_bf = (q_pp + q_p - cross(q_pp, q_p))/(1 - dot(q_pp, q_p));
% crp = [-0.443396226415094, 0, -0.103773584905660]


%% CC17
% 3
clear; clc;
sigma = [0.1, 0.2, 0.3];
sigma_sq = norm(sigma)^2;
sigma_shadow = -sigma./sigma_sq;
% mrp = [-0.714285714285714, -1.42857142857143, -2.14285714285714]

%% CC18
% 1
clear; clc;
mrp = [0.1, 0.2, 0.3];
dcm = mrp2dcm(mrp);
function dcm = mrp2dcm(mrp)
    mrp_sq = norm(mrp)^2;
    mrp_tilde = skew_symmetric(mrp);
    dcm = eye(3) + (8 * mrp_tilde * mrp_tilde - 4 * (1-mrp_sq) * mrp_tilde)/(1+mrp_sq)^2;
end

% 2
clear; clc
dcm = [0.763314, 0.0946746, -0.639053;
        -0.568047,-0.372781, -0.733728;
        -0.307692, 0.923077, -0.230769];
mrp = dcm2mrp(dcm);
% mrp = [-0.5, 0.1, 0.2 ]

function mrp = dcm2mrp(dcm)
    zeta = sqrt(trace(dcm) + 1);
    mrp = 1/(zeta*(zeta+2)) * [dcm(2,3)-dcm(3,2); dcm(3,1)-dcm(1,3); dcm(1,2)-dcm(2,1)];
end

%% CC19
% 2
clear; clc;
sigma_BN = [0.1, 0.2, 0.3];
sigma_RB = [-0.1, 0.3, 0.1];
dcm_BN = mrp2dcm(sigma_BN);
dcm_RB = mrp2dcm(sigma_RB);
dcm_RN = dcm_RB * dcm_BN;
sigma_RN = dcm2mrp(dcm_RN);
% mrp = [-0.1602, 0.4162, 0.5296 ]

% 3
clear; clc;
sigma_BN = [0.1, 0.2, 0.3];
sigma_RN = [0.5, 0.3, 0.1];
dcm_BN = mrp2dcm(sigma_BN);
dcm_RN = mrp2dcm(sigma_RN);
dcm_BR = dcm_BN * dcm_RN';
sigma_BR = dcm2mrp(dcm_BR);
% mrp = [-0.38, 0.1144, -0.0233 ]


%% CC2 - TRIAD
% 1
clear; clc;
v1_B = [0.8273, 0.5541, -0.0920];
v2_B = [-0.8285, 0.5522, -0.0955];

v1_N = [-0.1517, -0.9669, 0.2050];
v2_N = [-0.8393, 0.4494, -0.3044];

t1_B = v1_B;
v1bxv2b = cross(v1_B, v2_B);
t2_B = (v1bxv2b)/(norm(v1bxv2b));
t3_B = cross(t1_B, t2_B);
BT = [t1_B', t2_B', t3_B'];

t1_N = v1_N;
v1nxv2n = cross(v1_N, v2_N);
t2_N = (v1nxv2n)/(norm(v1nxv2n));
t3_N = cross(t1_N, t2_N);
NT = [t1_N', t2_N', t3_N'];

BN = BT * NT';
% dcm = [0.415527462219241, -0.855026538049592, 0.310026053309796,
%        -0.833866918551016, -0.494241815751786, -0.245448150461420,
% 0.363124807552830, -0.156564868983612, -0.918489856168244]

% 2
clear; clc

BN_est = [0.969846, 0.17101, 0.173648;
    -0.200706, 0.96461, 0.17101;
    -0.138258, -0.200706, 0.969846];
BN = [0.963592, 0.187303, 0.190809;
    -0.223042 0.956645, 0.187303;
    -0.147454, -0.223042, 0.963592];

C = BN_est * BN';
phi_C = acosd(0.5 * (C(1,1) + C(2,2) + C(3,3) - 1));
% C = BN_est * BN;

phi_BN = acosd(0.5 * (BN(1,1) + BN(2,2) + BN(3,3) - 1));
phi_BN_est = acosd(0.5 * (BN_est(1,1) + BN_est(2,2) + BN_est(3,3) - 1));
error = phi_BN - phi_BN_est;