function state_dot = TBP_duration(time, state, mu, T, I_sp, g_0, angle, delay_time, duration)
    % 2 Body Problem EOMs
    x = state(1);
    y = state(2);
    z = state(3);
    xdot = state(4);
    ydot = state(5);
    zdot = state(6);
    mass = state(7);

    r = sqrt(x^2 + y^2 + z^2);
    v = sqrt(xdot^2 + ydot^2 + zdot^2);

    v1 = [xdot, ydot, zdot]/v;
    v2 = cross([x, y, z]/r, v1);
    v3 = cross(v1, v2);

    % DCM_VNB_N
    VN = [v1; v2; v3];

    % Thrust only in the orbit normal 
    % u_thrust_VNB = [0, 1, 0]';
    u_thrust_VNB = [sin(angle), cos(angle), 0]';
    u_thrust_N = VN' * u_thrust_VNB;

    state_dot(1, 1) = xdot;
    state_dot(2, 1) = ydot;
    state_dot(3, 1) = zdot;
    
    if (time > delay_time && time < (delay_time + duration))
        state_dot(4, 1) = -mu/r^3*x + T/mass * u_thrust_N(1);
        state_dot(5, 1) = -mu/r^3*y + T/mass * u_thrust_N(2);
        state_dot(6, 1) = -mu/r^3*z + T/mass * u_thrust_N(3);
    
        state_dot(7, 1) = -T/(I_sp * g_0);
    else
        state_dot(4, 1) = -mu/r^3*x;
        state_dot(5, 1) = -mu/r^3*y;
        state_dot(6, 1) = -mu/r^3*z;
    
        state_dot(7, 1) = 0;
    end
end