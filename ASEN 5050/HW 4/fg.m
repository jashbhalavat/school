function rv = fg(r0, v0, theta_star_1, mu)
    % Calculate rv at a different position giving initial r0, v0, final theta_star, and
    % gravitational parameter of central body in a 2 BP

    coe = rv2coe(r0, v0, mu);
    h = norm(cross(r0, v0));
    p = h^2/mu;
    e = coe(2);
    theta_star_0 = coe(6);
    
    delta_theta_star = theta_star_1 - theta_star_0;

    r = p / (1 + e*cos(theta_star_1));

    r0_mag = norm(r0);

    f = 1 - r/p * (1 - cos(delta_theta_star));
    g = r0_mag*r/sqrt(mu*p) * sin(delta_theta_star);

    f_dot = sqrt(mu/p) * tan(delta_theta_star/2) * ((1 - cos(delta_theta_star))/(p) - 1/r - 1/r0_mag);
    g_dot = 1 - (r0_mag/p) * (1 - cos(delta_theta_star));

    r_out = f*r0 + g*v0;
    v_out = f_dot*r0 + g_dot*v0;

    rv = [r_out; v_out];
end

