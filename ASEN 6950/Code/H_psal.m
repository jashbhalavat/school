function out = H_psal(V, V_star, delta_s, n_hat, statef)
    out = [statef(1) - V(1);
            statef(2) - V(2);
            statef(3) - V(3);
            statef(4) - V(4);
            statef(6) - V(6);
            V(2);
            dot(V-V_star, n_hat) - delta_s];
end