function state_phi_dot = CR3BP_full(state_phi, mu, F, Isp, u_hat, t_star, l_star, M_sc_0)
    % Full state vector and state transition matrix differential equation
    % Inputs:
    % state_phi - Augmented state vector and STM [42x1]. The state vector -
    % [x0, y0, z0, x0_dot, y0_dot, z0_dot, mass]. The STM - is 7x7 with each
    % element described as - phi_ij = dxi(tf)/dxj(t0). The phi matrix is
    % reshaped such that all the rows are concatenated vertically. For
    % example - 
    % phi_mat = [phi11, phi12, phi13, ..., phi17;
    %           [phi21, phi22, phi23, ..., phi27;
    %           ...
    %           [phi71, phi72, phi73, ..., phi77]
    % becomes
    % phi_row = [phi11, phi12, ..., phi17, phi21, phi22, ..., phi77]'
    % 
    % mu - system mass ratio [-]
    % 
    % Output
    % state_phi_dot - Augmented state vector dot and STM_dot [42x1]. The
    % augmentation and reshaping scheme remains the same as the input.

    x = state_phi(1);
    y = state_phi(2);
    z = state_phi(3);
    xdot = state_phi(4);
    ydot = state_phi(5);
    zdot = state_phi(6);
    mass = state_phi(7);

    r1 = sqrt((x + mu)^2 + (y)^2 + (z)^2);
    r2 = sqrt((x - 1 + mu)^2 + (y)^2 + (z)^2);

    % Non-dim thrust
    f = (F*t_star^2)/(l_star*M_sc_0);
    % Non-dim mass
    m = M_sc/M_sc_0;
    % acceleration due to low thrust
    a_lt = f/m*u_hat;

    state_dot(1, 1) = xdot;
    state_dot(2, 1) = ydot;
    state_dot(3, 1) = zdot;

    state_dot(4, 1) = 2*ydot + x - (1 - mu)*(x + mu)/(r1^3) - mu * (x - 1 + mu)/(r2^3) + a_lt(1);
    state_dot(5, 1) = -2*xdot + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3) + a_lt(2);
    state_dot(6, 1) = -(1 - mu)*z/(r1^3) - mu*z/(r2^3) + a_lt(3);

    % g is assumed to be 9.80665e-3 km/s^2
    mdot = -(f*l_star)/(Isp*9.80665e-3*t_star);
    state_dot(7, 1) = mdot;
    
    % Calc pseudo-potentials
    uxx = u_xx(mu, [x, y, z]);
    uyy = u_yy(mu, [x, y, z]);
    uxy = u_xy(mu, [x, y, z]);
    uzz = u_zz(mu, [x, y, z]);
    uxz = u_xz(mu, [x, y, z]);
    uyz = u_yz(mu, [x, y, z]);

    U_mat = [uxx, uxy uxz; uxy, uyy uyz; uxz uyz uzz];
    Omega = [0 2 0; -2 0 0; 0 0 0];
    A = [zeros(3), eye(3);
        U_mat, Omega];

    % Get only the phi elements into a row
    phi_row = state_phi(8:end);

    % Converting phi to matrix
    phi_mat = reshape(phi_row, [7,7])';

    % Get phi_dot
    phi_dot_mat = A * phi_mat;

    % Convert back to row
    phi_dot_row = reshape(phi_dot_mat', [42,1]);

    % Augment state and phi (in row form)
    state_phi_dot = [state_dot; phi_dot_row];

end