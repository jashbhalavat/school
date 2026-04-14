function DF = DF_mat(V, statef_V, options, system_params)
    % Modified constraint DF matrix

    % Get mass ratio of system
    mu = system_params(1);
    t_star = system_params(2);
    l_star = system_params(3);
    T = system_params(4);
    Isp = system_params(5);
    M_sc_0 = system_params(6);
    f = system_params(7);
    mdot = system_params(8);

    x_1_0 = V(1:6);
    m_1 = V(7);
    Dt1 = V(8);
    x_2_0 = V(9:14);
    m_2_0 = V(15);
    uhat_2 = V(16:18);
    Dt2 = V(19);
    x_3_0 = V(20:25);
    m_3_0 = V(26);
    uhat_3 = V(27:29);
    Dt3 = V(30);
    x_4_0 = V(31:36);
    m_4 = V(37);
    Dt4 = V(38);

    x_1_f = statef_V(1:6);
    x_2_f = statef_V(9:14);
    m_2_f = statef_V(15);
    x_3_f = statef_V(20:25);
    m_3_f = statef_V(26);
    x_4_f = statef_V(31:36);

    phi0_36 = reshape(eye(6), [36, 1]); % Initial phi is identity
    phi0_100 = reshape(eye(10), [100, 1]); % Initial phi is identity

    x1_0 = [V(1:6); phi0_36];
    x2_0 = [V(9:18); phi0_100];
    x3_0 = [V(20:29); phi0_100];
    x4_0 = [V(31:36); phi0_36];

    % Define functions
    init_fun = @(t,state)CR3BP_full_mass(state, mu, T, Isp, t_star, l_star, M_sc_0);
    final_fun = @(t,state)CR3BP_full_mass(state, mu, T, Isp, t_star, l_star, M_sc_0);

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
    xdot_2_f = CR3BP_with_non_dim_mass(x2_f, mu, uhat_2, f, mdot);
    mdot_2_f = xdot_2_f(7);
    xdot_2_f = xdot_2_f(1:6);
    xdot_3_f = CR3BP_with_non_dim_mass(x3_f, mu, uhat_3, f, mdot);
    mdot_3_f = xdot_3_f(7);
    xdot_3_f = xdot_3_f(1:6);
    xdot_4_0 = CR3BP(x_4_0, mu);
    xdot_4_f = CR3BP(x4_f, mu);

    phi_row_1 = x1_f(7:end);
    phi_mat_1 = reshape(phi_row_1, [6,6])';
    phi_row_2 = x2_f(11:end);
    phi_mat_2 = reshape(phi_row_2, [10,10])';
    phi_row_3 = x3_f(11:end);
    phi_mat_3 = reshape(phi_row_3, [10,10])';
    phi_row_4 = x4_f(7:end);
    phi_mat_4 = reshape(phi_row_4, [6,6])';

    % out = [phi_mat_1, xdot_1_f, -eye(6), zeros(6,5), zeros(6,11), zeros(6,8);
    %         zeros(9,7), [phi_mat_2(1:6, 1:6); zeros([2,6]); phi_mat_2(7, 1:6)], [phi_mat_2(1:6, 7); 0; 1; phi_mat_2(7,7)], [xdot_2_f; 0; 0; mdot_2_f], [phi_mat_2(1:6, 8:10); [2*uhat_2(1), 2*uhat_2(2), 2*uhat_2(3)]; zeros(1,3); phi_mat_2(7, 8:10)], [-eye(6), zeros(6,5); zeros(2,11); zeros(1,6), -1, zeros(1,4)], zeros(9,8);
    %         zeros(8,18), [phi_mat_3(1:6, 1:7); zeros(1,7); zeros(1,6), 1], [xdot_3_f; 0; mdot_3_f], [phi_mat_3(1:6, 8:10); [2*uhat_3(1), 2*uhat_3(2), 2*uhat_3(3)]; zeros(1,3)], [-eye(6); zeros(2,6)], [zeros(6,1); 0; -1], [zeros(8,1)];
    %         zeros(6,29), [phi_mat_4, zeros(6,1), xdot_4_f]];

    % F1V1 = [eye(6), zeros(6,2);
    %         phi_mat_1, zeros(6,1), xdot_1_f];
    % F1V2 = [zeros(6,11);
    %         -eye(6), zeros(6,5)];
    % F1V3 = zeros(12,11);
    % F1V4 = zeros(12,8);

    F1V1 = [eye(7), zeros(7,1);
                phi_mat_1, zeros(6,1), xdot_1_f
                zeros(1,6), 1, 0];
    F1V2 = [zeros(7,11);
                -eye(6), zeros(6,5)
                zeros(1,6), -1, zeros(1,4)];
    F1V3 = zeros(14,11);
    F1V4 = zeros(14,8);

    F2V1 = zeros(8);
    F2V2 = [phi_mat_2(1:7,:), [xdot_2_f; mdot_2_f];
            zeros(1,7), 2*uhat_2', 0];
    F2V3 = [-eye(6), zeros(6,5);
            zeros(1,6), -1, zeros(1,4);
            zeros(1,11)];
    F2V4 = zeros(8,8);

    F3V1 = zeros(8);
    F3V2 = zeros(8,11);
    F3V3 = [phi_mat_3(1:7,:), [xdot_3_f; mdot_3_f];
            zeros(1,7), 2*uhat_3', 0];
    F3V4 = [-eye(6), zeros(6,2);
            zeros(2,8)];
    F3V4(7,7) = -1;

    F4V1 = zeros(6,8);
    F4V2 = zeros(6,11);
    F4V3 = zeros(6,11);
    F4V4 = [phi_mat_4, zeros(6,1), xdot_4_f];

    DF = [F1V1, F1V2, F1V3, F1V4;
        F2V1, F2V2, F2V3, F2V4;
        F3V1, F3V2, F3V3, F3V4;
        F4V1, F4V2, F4V3, F4V4];

end