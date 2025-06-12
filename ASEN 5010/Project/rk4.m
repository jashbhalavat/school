function [x, u] = rk4(fn, t, x0, I_B, u)

    dt = t(2) - t(1);
    
    % x is nx6
    x = x0';

    for i = 1:length(t)
        k1 = dt * fn(x(i,:), t(i), u, I_B)';
        k2 = dt * fn(x(i,:) + k1/2, t(i) + dt/2, u, I_B)';
        k3 = dt * fn(x(i,:) + k2/2, t(i) + dt/2, u, I_B)';
        k4 = dt * fn(x(i,:) + k3, t(i) + dt, u, I_B)';

        x(i+1,:) = x(i,:) + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

        sigma_BN = x(i+1,1:3);
        if norm(sigma_BN) > 1
            x(i+1,1:3) = mrp_shadow(sigma_BN);
        end
    end
end



