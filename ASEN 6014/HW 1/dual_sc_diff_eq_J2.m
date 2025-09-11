function state_dot = dual_sc_diff_eq_J2(t, state)
    mu = 398600; % km3/s2
    J2 = 1082.63e-6;
    req = 6378.0; % km

    r1_N = state(1:3);
    v1_N = state(4:6);
    r2_N = state(7:9);
    v2_N = state(10:12);

    r1 = norm(r1_N);
    r2 = norm(r2_N);
    aj2_1 = -3/2 * J2 * mu/r1^2 * (req/r1)^2 * (1 - 5*(r1_N(3)/r1)^2) * r1_N/r1;
    aj2_2 = -3/2 * J2 * mu/r2^2 * (req/r2)^2 * (1 - 5*(r2_N(3)/r2)^2) * r2_N/r2;

    r1_dotdot = -mu/r1^3 * r1_N + aj2_1;
    r2_dotdot = -mu/r2^3 * r2_N + aj2_2;

    state_dot(1:3,1) = v1_N;
    state_dot(4:6,1) = r1_dotdot;
    state_dot(7:9,1) = v2_N;
    state_dot(10:12,1) = r2_dotdot;
end