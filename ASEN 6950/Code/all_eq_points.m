function [x_Ls, validity] = all_eq_points(mu)
    % Function that finds all equilibrium point positions
    
    % Initial guess for L1 is midpoint of P1 and P2
    x_L1_0 = (-mu + (1 - mu))/2;
    
    % Call function to get position of L1 and validity of L1
    [x_L1, x_L1_validity] = find_coll_eq_pts(x_L1_0, mu);
    
    % Initial guess for L2 is distance between P2 and L1. That is added to P2
    % position.
    x_L2_0 = (1-mu)-x_L1 + (1-mu);
    
    % Call function to get position of L2 and validity of L2
    [x_L2, x_L2_validity] = find_coll_eq_pts(x_L2_0, mu);
    
    % Initial guess for L3 is distance between P2 and L1. That is subtracted
    % from P1 position.
    x_L3_0 = -mu - ((1-mu)-x_L1);
    
    % Call function to get position of L2 and validity of L2
    [x_L3, x_L3_validity] = find_coll_eq_pts(x_L3_0, mu);

    % Find L4 and L5 points
    x_L4 = 1/2 - mu;
    y_L4 = sqrt(3)/2;
    y_L5 = -sqrt(3)/2;

    x_Ls = [x_L1, 0; x_L2, 0; x_L3, 0; x_L4, y_L4; x_L4, y_L5];
    validity = [x_L1_validity, x_L2_validity, x_L3_validity, true, true];
end
