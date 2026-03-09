function out = DF_mat(V, statef_V, options, system_params, x_1_des, x_4_des)
    % Modified constraint DF matrix

    % Get mass ratio of system
    mu = system_params(1);
    t_star = system_params(2);
    l_star = system_params(3);
    T = system_params(4);
    Isp = system_params(5);
    M_sc_0 = system_params(6);

    x_1_0 = V(1:6);
    Dt1 = V(7);
    x_2_0 = V(8:13);
    m_2_0 = V(14);
    Dt2 = V(15);
    uhat_2 = V(16:18);
    x_3_0 = V(19:24);
    m_3_0 = V(25);
    Dt3 = V(26);
    uhat_3 = V(27:29);
    x_4_0 = V(30:35);
    Dt4 = V(36);

    x_1_f = statef_V(1:6);
    x_2_f = statef_V(8:13);
    m_2_f = statef_V(14);
    x_3_f = statef_V(19:24);
    m_3_f = statef_V(25);
    x_4_f = statef_V(30:35);

    phi0 = reshape(eye(6), [36, 1]); % Initial phi is identity

    x1_0 = [V(1:6); phi0];
    x2_0 = [V(8:13); m_2_0; phi0];
    x3_0 = [V(15:20); m_3_0; phi0];
    x4_0 = [V(22:27); phi0];

    % Define functions
    init_fun = @(t,state)CR3BP_full_mass(state, mu, T, Isp, uhat_2, t_star, l_star, m_2_0);
    final_fun = @(t,state)CR3BP_full_mass(state, mu, T, Isp, uhat_3, t_star, l_star, m_3_0);

    % Non-dim thrust
    f = (T*t_star^2)/(l_star*M_sc_0);
    % Non-dim mass
    m_2 = m_2_0/M_sc_0;
    % acceleration due to low thrust
    a_lt_2 = -f/m_2*uhat_2;
    a_lt_mat_2 = [a_lt_2(1), 0, 0; 0, a_lt_2(2), 0; 0, 0, a_lt_2(3)];
    a_lt_dot_2 = -f/m_2^2*uhat_2;

    m_3 = m_3_0/M_sc_0;
    a_lt_3 = -f/m_3*uhat_3;
    a_lt_mat_3 = [a_lt_3(1), 0, 0; 0, a_lt_3(2), 0; 0, 0, a_lt_3(3)];
    a_lt_dot_3 = -f/m_3^2*uhat_3;

    [~, x1_out] = ode113(@(t,state)CR3BP_full(state, mu), [0, Dt1], x1_0, options);
    [~, x2_out] = ode113(init_fun, [0 Dt2], x2_0, options);
    [~, x3_out] = ode113(final_fun, [0 Dt3], x3_0, options);
    [~, x4_out] = ode113(@(t,state)CR3BP_full(state, mu), [0, Dt4], x4_0, options);
    
    x1_f = x1_out(end, :)';
    x2_f = x2_out(end, :)';
    x3_f = x3_out(end, :)';
    x4_f = x4_out(end, :)';
    xdot_1_0 = CR3BP(x_1_0, mu);
    xdot_1_f = CR3BP(x1_f, mu);
    xdot_2_f = CR3BP_with_non_dim_mass(x2_f, mu, T, Isp, uhat_2, t_star, l_star, m_2_0);
    mdot_2_f = xdot_2_f(7);
    xdot_2_f = xdot_2_f(1:6);
    xdot_3_f = CR3BP_with_non_dim_mass(x3_f, mu, T, Isp, uhat_3, t_star, l_star, m_3_0);
    mdot_3_f = xdot_3_f(7);
    xdot_3_f = xdot_3_f(1:6);
    xdot_4_f = CR3BP(x4_f, mu);
    xdot_1_des = CR3BP(x_1_des, mu);
    xdot_4_des = CR3BP(x_4_des, mu);

    phi_row_1 = x1_f(7:end);
    phi_mat_1 = reshape(phi_row_1, [6,6])';
    phi_row_2 = x2_f(8:end);
    phi_mat_2 = reshape(phi_row_2, [6,6])';
    phi_row_3 = x3_f(8:end);
    phi_mat_3 = reshape(phi_row_3, [6,6])';
    phi_row_4 = x4_f(7:end);
    phi_mat_4 = reshape(phi_row_4, [6,6])';

    out = [phi_mat_1, xdot_1_f, -eye(6), zeros([6,23]);
            zeros(8,7), [phi_mat_2; zeros([2,6])], [zeros([3,1]); a_lt_dot_2; zeros(2,1)], [xdot_2_f; mdot_2_f; 0], [zeros(3); a_lt_mat_2; zeros(1,3); [2*uhat_2(1), 2*uhat_2(2), 2*uhat_2(3)]], [-eye(6); zeros(2,6)], [zeros(6,1); -1; 0], zeros([8,11]);
            zeros(7,18), [phi_mat_3; zeros(1,6)], [zeros(3,1); a_lt_dot_3; 0], [xdot_3_f; 0], [zeros(3); a_lt_mat_3; [2*uhat_3(1), 2*uhat_3(2), 2*uhat_3(3)]], [-eye(6); zeros(1,6)], zeros(7,1);
            [eye(6); zeros(6)], [xdot_1_0 - xdot_1_des; zeros(6,1)], zeros(12,22), [zeros(6); eye(6)], [zeros(6,1); xdot_4_f - xdot_4_des]];

end