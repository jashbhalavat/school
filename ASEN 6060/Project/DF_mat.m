function out = DF_mat(V, options, mu)
    % Modified constraint DF matrix

    phi0 = reshape(eye(6), [36, 1]); % Initial phi is identity

    x1_0 = [V(1:6); phi0];
    x2_0 = [V(8:13); phi0];
    x3_0 = [V(15:20); phi0];
    x4_0 = [V(22:27); phi0];

    [~, x1_out] = ode113(@(t,state)CR3BP_full(state, mu), [0 V(7)], x1_0, options);
    [~, x2_out] = ode113(@(t,state)CR3BP_full(state, mu), [0 V(14)], x2_0, options);
    [~, x3_out] = ode113(@(t,state)CR3BP_full(state, mu), [0 V(21)], x3_0, options);
    [~, x4_out] = ode113(@(t,state)CR3BP_full(state, mu), [0 V(28)], x4_0, options);

    x1_f = x1_out(end, :)';
    x2_f = x2_out(end, :)';
    x3_f = x3_out(end, :)';
    x4_f = x4_out(end, :)';

    phi_row_1 = x1_f(7:end);
    phi_mat_1 = reshape(phi_row_1, [6,6])';
    phi_row_2 = x2_f(7:end);
    phi_mat_2 = reshape(phi_row_2, [6,6])';
    phi_row_3 = x3_f(7:end);
    phi_mat_3 = reshape(phi_row_3, [6,6])';
    phi_row_4 = x4_f(7:end);
    phi_mat_4 = reshape(phi_row_4, [6,6])';

    out = [phi_mat_1(1:3,:), x1_f(4:6), -eye([3,6]), zeros([3,15]);
            zeros([3,7]), phi_mat_2(1:3,:), x2_f(4:6), -eye([3,6]), zeros([3,8]);
            zeros([3,14]), phi_mat_3(1:3,:), x3_f(4:6), -eye([3,6]), zeros([3,1]);
            zeros([3,21]), phi_mat_4(1:3,:), x4_f(4:6)];

end