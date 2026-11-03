function fx = collinear_eq_pts(x, mu)
    % Function to get roots for collinear equilibrium points
    fx = x - ((1-mu)*(x+mu))/(abs(x+mu)^3) - (mu*(x-1+mu))/(abs(x-1+mu)^3);
end