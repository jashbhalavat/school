function coe = rv2coe(r, v, mu)
    % Convert pos/vel to classical orbital elements
    % Output vector is arranged as follows:
    % coe - [a, e, i, raan, aop, theta_star]
    
    r_norm = norm(r);
    r_hat = r/r_norm;
    
    v_norm = norm(v);
    v_hat = v/v_norm;
    
    h = cross(r, v);
    h_norm = norm(h);
    h_hat = h/h_norm;
    
    K = [0 0 1];
    n = cross(K, h);
    n_norm = norm(n);
    n_hat = n/n_norm;
    
    e = 1/mu * ((v_norm^2 - mu/r_norm)*r - dot(r,v)*v);
    e_norm = norm(e);
    e_hat = e/e_norm;
    
    eps = v_norm^2/2 - mu/r_norm;
    
    if e == 1.0
        p = h_norm^2/mu;
    else
        a = -mu/(2*eps);
        p = a * (1 - e_norm^2);
    end
    
    i = acos(h_hat(3));
    raan = sign(n_hat(2)) * abs(acos(n_hat(1)));
    aop = sign(e(3)) * abs(acos(dot(n, e)/(n_norm * e_norm)));
    theta_star = sign(dot(r, v)) * abs(acos(dot(e, r)/(e_norm * r_norm)));
    
    % Special cases
    aop_true = acos(e_hat(1));
    if e(2) < 0
        aop_true = 2*pi - aop_true;
    end

    coe = [a, e_norm, i, raan, aop, theta_star];
end