function out = DF_mat(V, options, mu, T, Isp, uncorrected_init_thrust, uncorrected_final_thrust, t_star_em, l_star_em, init_mass, uncorrected_final_mass, desired_stage_1, desired_stage_4)
    % Modified constraint DF matrix

    phi0 = reshape(eye(6), [36, 1]); % Initial phi is identity

    xl1_0 = [V(1:6); phi0];
    x1_0 = [V(8:13); init_mass; phi0];
    x2_0 = [V(15:20); uncorrected_final_mass; phi0];
    xl2_0 = [V(22:27); phi0];

    init_fun = @(t,state)CR3BP_full_mass(state, mu, T, Isp, uncorrected_init_thrust, t_star_em, l_star_em, init_mass);
    final_fun = @(t,state)CR3BP_full_mass(state, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, init_mass);

    [~, xl1_out] = ode113(@(t,state)CR3BP_full(state, mu), [0, V(7)], xl1_0, options);
    [~, x1_out] = ode113(init_fun, [0 V(14)], x1_0, options);
    [~, x2_out] = ode113(final_fun, [0 V(21)], x2_0, options);
    [~, xl2_out] = ode113(@(t,state)CR3BP_full(state, mu), [0, V(28)], xl2_0, options);
    
    xl1_0 = xl1_out(1, :)';
    xl1_f = xl1_out(end, :)';
    x1_f = x1_out(end, :)';
    x2_f = x2_out(end, :)';
    xl2_f = xl2_out(end, :)';
    xdot_l1_0 = CR3BP(xl1_0, mu);
    xdot_l1_f = CR3BP(xl1_f, mu);
    xdot_1_f = CR3BP_with_non_dim_mass(x1_f, mu, T, Isp, uncorrected_init_thrust, t_star_em, l_star_em, init_mass);
    xdot_2_f = CR3BP_with_non_dim_mass(x2_f, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, init_mass);
    xdot_l2_f = CR3BP(xl2_f, mu);
    desired_stage_1_dot = CR3BP(desired_stage_1, mu);
    desired_stage_4_dot = CR3BP(desired_stage_4, mu);

    phi_row_l1 = xl1_f(7:end);
    phi_mat_l1 = reshape(phi_row_l1, [6,6])';
    phi_row_1 = x1_f(8:end);
    phi_mat_1 = reshape(phi_row_1, [6,6])';
    phi_row_2 = x2_f(8:end);
    phi_mat_2 = reshape(phi_row_2, [6,6])';
    phi_row_l2 = xl2_f(7:end);
    phi_mat_l2 = reshape(phi_row_l2, [6,6])';


    % out = [phi_mat_1, xdot_1_f(1:6), -eye(6), zeros([6,1])];

    out = [phi_mat_l1, xdot_l1_f, -eye(6), zeros([6,15]);
        zeros([6,7]), phi_mat_1, xdot_1_f(1:6), -eye(6), zeros([6,8]);
        zeros([6,14]), phi_mat_2, xdot_2_f(1:6), -eye(6), zeros([6,1]);
        zeros([6,21]), phi_mat_l2, xdot_l2_f - desired_stage_4_dot;
        eye(6), xdot_l1_0 - desired_stage_1_dot, zeros([6,21])];

end