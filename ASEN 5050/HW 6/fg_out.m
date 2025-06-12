function out = fg_out(r0, r1, delta_theta_star, p, mu)
    % Calculate rv at a different position giving initial r0, v0, final theta_star, and
    % gravitational parameter of central body in a 2 BP

    r = norm(r1);

    r0_mag = norm(r0);

    f = 1 - r/p * (1 - cos(delta_theta_star));
    g = r0_mag*r/sqrt(mu*p) * sin(delta_theta_star);

    f_dot = sqrt(mu/p) * tan(delta_theta_star/2) * ((1 - cos(delta_theta_star))/(p) - 1/r - 1/r0_mag);
    g_dot = 1 - (r0_mag/p) * (1 - cos(delta_theta_star));

    out = [f g f_dot g_dot];
end

