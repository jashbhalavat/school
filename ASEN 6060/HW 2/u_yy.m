function out = u_yy(mu, x_eq)
    % Pseudo potential function partial derivative wrt y, y at eq point
    % Assuming z = 0
    x = x_eq(1);
    y = x_eq(2);
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2);
    out = 1 - (1-mu)/(r1^3) - mu/(r2^3) + (3*(1-mu)*y^2)/(r1^5) + (3*mu*y^2)/(r2^5);
end
