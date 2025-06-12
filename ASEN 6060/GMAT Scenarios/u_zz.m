function out = u_zz(mu, x_eq)
    % Pseudo potential function partial derivative wrt z, z at eq point
    % Assuming z = 0
    x = x_eq(1);
    y = x_eq(2);
    z = x_eq(3);
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    out = -(1-mu)/(r1^3) - mu/(r2^3) + (3 * (1-mu) * z^2)/(r1^5) + (3*mu*z^2)/(r2^5);
end