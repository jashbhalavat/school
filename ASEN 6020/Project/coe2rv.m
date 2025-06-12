function rv = coe2rv(coe, mu)
    % Compute pos/vel from classical orbital elements where coe is:
    % coe - [a, e, i, raan, aop, theta_star]

    a = coe(1);
    e = coe(2);
    i = coe(3);
    raan = coe(4);
    aop = coe(5);
    theta_star = coe(6);

    h = sqrt(a * (1 - e^2) * mu);

    v_r = mu/h * e * sin(theta_star);
    v_theta = mu/h * (1 + e*cos(theta_star));

    r = (h^2/mu)/(1 + e*cos(theta_star));
    r_rth = [r, 0, 0];

    v_rth = [v_r, v_theta, 0];

    theta = theta_star + aop;

    r_xyz = R3(-raan) * R1(-i) * R3(-theta) * r_rth';
    v_xyz = R3(-raan) * R1(-i) * R3(-theta) * v_rth';

    rv = [r_xyz; v_xyz];
end


