function [x_Lx, validity] = find_coll_eq_pts(x0, mu)
    % Function to find collinear equilibrium points

    % Function for collinear equilibrium points
    fun = @(x)collinear_eq_pts(x, mu);

    % Set display to zero and thresholds to 1e-17 to get as close to double
    % precision as possible
    options = optimoptions('fsolve', 'Display','none', 'FunctionTolerance',1e-17, ...
        'StepTolerance',1e-17, 'OptimalityTolerance',1e-17);
    
    % Use fsolve to find root of function
    x_Lx = fsolve(fun, x0, options);

    % Test if the root is sufficiently close to 0. If so, validity is true,
    % else validity is false
    test_Lx = collinear_eq_pts(x_Lx, mu);
    validity = false;
    if abs(test_Lx) < 1e-15
        validity = true;
    end
end