function V_family = natural_param_continuation(V_soln, mu)
    % V_soln - First corrected free variable solution
    % mu - system parameter

    % out - Family of free variables
    % ASSUMPTION - Earth-Moon System!

    a = 384400; % [km] Average earth-moon semi-major axis
    r_E = 6378.1363; % [km] Earth equatorial radius (Vallado, Appendix D)
    r_Moon = 1738; % [km] Moon equatorial radius (Vallado, Appendix D)
    r_E_normalized = r_E/a;
    r_Moon_normalized = r_Moon/a;

    TOL = 5e-14;

    V_family = V_soln;
    delta = -2e-3;

    % Set options for ode113
    options = odeset('RelTol', TOL, 'AbsTol', TOL);

    [tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

    % Check if initial position is outside the P1, P2 bodies' normalized radius
    initial_x_pos = V_family(1:3,1);
    p1_pos = [-mu, 0, 0]';
    p2_pos = [1-mu, 0, 0]';
    p1_minus_init_pos = p1_pos - initial_x_pos;
    p2_minus_init_pos = p2_pos - initial_x_pos;
    p1_or_p2_bound = false;
    if (norm(p1_minus_init_pos) <= r_E_normalized)
        p1_or_p2_bound = true;
    end
    if (norm(p2_minus_init_pos) <= r_E_normalized)
        p1_or_p2_bound = true;
    end

    family_member = 1;
    max_family_members = 150;

    while ((~p1_or_p2_bound) && (family_member <= max_family_members))
        print_out = sprintf('Family member - %d', family_member);
        disp(print_out)
        V = V_family(:,family_member);
        x_star = V(1);
        x_star_plus_delta = x_star + delta;
        H_norm = norm(H(V_soln(1:6), xout_corrected(end,:), x_star_plus_delta));
        counter = 1;

        % Check if initial position is outside the P1, P2 bodies' normalized radius
        initial_x_pos = V_family(1:3,family_member);
        p1_minus_init_pos = p1_pos - initial_x_pos;
        p2_minus_init_pos = p2_pos - initial_x_pos;
        p1_or_p2_bound = false;
        if (norm(p1_minus_init_pos) <= r_E_normalized)
            p1_or_p2_bound = true;
        end
        if (norm(p2_minus_init_pos) <= r_E_normalized)
            p1_or_p2_bound = true;
        end
        
        while ((H_norm > TOL) && (counter < 10))
            V_i = V(:,counter);
            phi_0 = reshape(eye(6), [36,1]);
            state0 = [V_i(1:6); phi_0];

            [tout, state_out] = ode113(@(t, state)CR3BP_full(state, mu), [0, V_i(end)], state0, options);
            
            statef = state_out(end, :)';
            state_dot = CR3BP_full(statef, mu);
            H_V = H(V_i(1:6), statef, x_star_plus_delta);
            H_norm = norm(H_V);
            DH_V = DH_mat(statef, state_dot);
            try
                inv(DH_V);
            catch
                break
            end
            V_ip1 = V_i - inv(DH_V) * H_V;
            V(:,counter+1) = V_ip1;
            counter = counter + 1;
        end

        V_family(:,family_member+1) = V(:,end);
        family_member = family_member + 1;
    end
end
    