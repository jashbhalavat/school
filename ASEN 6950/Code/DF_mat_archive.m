function out = DF_mat(state, state_dot)
    % Modified constraint DF matrix

    % Convert STM to matrix from row vector
    phi_row = state(7:end);
    phi_mat = reshape(phi_row, [6,6])';
    % Subtract identity and extract all rows except 5th row
    phi_minus_I = phi_mat - eye(6);
    phi_mod = [phi_minus_I(1:4, :); phi_minus_I(6,:)];

    % Grab states - x_dot, ydot, zdot, xdotdot, zdotdot
    state_dot_col = [state_dot(1:4); state_dot(6)];

    out = [phi_mod, state_dot_col; 0 1 0 0 0 0 0];
end