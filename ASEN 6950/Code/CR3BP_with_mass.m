function state_dot = CR3BP_with_mass(state, mu, T, Isp, T_star, u_rot)
    % Circular Restricted 3 Body Problem non-dimensional EOMs
    % Inputs
    %   state - [x, y, z, xdot, ydot, zdot, mass]
    %   mu - earth-moon non-dim mass ratio
    %   T - thrust [Nm]
    %   Isp - Specific impulse [sec]
    % Output
    %   state_dot - [xdot, ydot, zdot, xdotdot, ydotdot, zdotdot, mass_dot]

    x = state(1);
    y = state(2);
    z = state(3);
    xdot = state(4);
    ydot = state(5);
    zdot = state(6);

    r1 = sqrt((x + mu)^2 + (y)^2 + (z)^2);
    r2 = sqrt((x - 1 + mu)^2 + (y)^2 + (z)^2);

    state_dot(1, 1) = xdot;
    state_dot(2, 1) = ydot;
    state_dot(3, 1) = zdot;

    % S/C velocity unit vector in rotating frame
    % v_rot = [state(4), state(5), state(6)];
    % v_hat_rot = v_rot/norm(v_rot);

    % Only firing in the velocity direction
    % u_rot = v_hat_rot;

    state_dot(4, 1) = 2*ydot + x - (1 - mu)*(x + mu)/(r1^3) - mu * (x - 1 + mu)/(r2^3) + T_star/state(7)*u_rot(1);
    state_dot(5, 1) = -2*xdot + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3) + T_star/state(7)*u_rot(2);
    state_dot(6, 1) = - (1 - mu)*z/(r1^3) - mu*z/(r2^3) + T_star/state(7)*u_rot(3);

    % g is assumed to be 9.81 m/s^2
    mdot = -T/(Isp*9.81);
    state_dot(7, 1) = mdot;
end