function out = DH_mat(state, state_dot)
    phi_row = state(7:end);
    phi_mat = reshape(phi_row, [6,6])';
    phi_minus_I = phi_mat - eye(6);
    phi_mod = [phi_minus_I(1:4, :); phi_minus_I(6,:)];

    state_dot_col = [state_dot(1:4); state_dot(6)];

    out = [phi_mod, state_dot_col; 0 1 0 0 0 0 0; 1 0 0 0 0 0 0];
end