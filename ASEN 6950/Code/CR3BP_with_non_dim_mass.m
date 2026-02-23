function state_dot = CR3BP_with_non_dim_mass(state, mu, F, Isp, u_hat, t_star, l_star, M_sc_0)
    % Circular Restricted 3 Body Problem non-dimensional EOMs
    % Source - https://ntrs.nasa.gov/api/citations/20180006443/downloads/20180006443.pdf
    % Inputs
    %   state - [x, y, z, xdot, ydot, zdot, mass] - Units [non-dim pos/vel, mass - kg]
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
    M_sc = state(7);

    r1 = sqrt((x + mu)^2 + (y)^2 + (z)^2);
    r2 = sqrt((x - 1 + mu)^2 + (y)^2 + (z)^2);

    state_dot(1, 1) = xdot;
    state_dot(2, 1) = ydot;
    state_dot(3, 1) = zdot;

    % Non-dim thrust
    f = (F*t_star^2)/(l_star*M_sc_0);
    % Non-dim mass
    m = M_sc/M_sc_0;
    % acceleration due to low thrust
    a_lt = f/m*u_hat;

    state_dot(4, 1) = 2*ydot + x - (1 - mu)*(x + mu)/(r1^3) - mu * (x - 1 + mu)/(r2^3) + a_lt(1);
    state_dot(5, 1) = -2*xdot + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3) + a_lt(2);
    state_dot(6, 1) = - (1 - mu)*z/(r1^3) - mu*z/(r2^3) + a_lt(3);

    % g is assumed to be 9.80665e-3 km/s^2
    mdot = -(f*l_star)/(Isp*9.80665e-3*t_star);
    state_dot(7, 1) = mdot;
end