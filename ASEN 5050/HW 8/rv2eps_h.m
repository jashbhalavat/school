function [eps, h_norm] = rv2eps_h(r, v, mu)
    % Convert pos/vel to eps and h
    % Output vector is arranged as follows:
    % [eps, h]
    
    r_norm = norm(r);
    r_hat = r/r_norm;
    
    v_norm = norm(v);
    v_hat = v/v_norm;
    
    h = cross(r, v);
    h_norm = norm(h);
    
    eps = v_norm^2/2 - mu/r_norm;
    
end