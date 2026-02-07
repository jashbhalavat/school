function out = u_yz(mu, x_eq)
    % Pseudo potential function partial derivative wrt y, z at eq point
    x = x_eq(1);
    y = x_eq(2);
    z = x_eq(3);
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    out = 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5;
end