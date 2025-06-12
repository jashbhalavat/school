function V_family = test_psal_cont(V_soln, mu)
    % V_soln - First corrected free variable solution
    % mu - system parameter

    % out - Family of free variables

    a = 384400; % [km] Average earth-moon semi-major axis
    r_E = 6378.1363; % [km] Earth equatorial radius (Vallado, Appendix D)
    r_Moon = 1738; % [km] Moon equatorial radius (Vallado, Appendix D)
    r_E_normalized = r_E/a;
    r_Moon_normalized = r_Moon/a;

    TOL = 1e-13;

    V_family = V_soln;
    % Arbitrarily pick an initial delta_s
    % delta_s = -0.01; % negative for a1
    % delta_s = -0.0005; % negative for a2
    % delta_s = 0.0005; % works for a2 pos
    % delta_s = -0.001; % works for a3
    delta_s = 1e-3;

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

    % Initialize nHat
    phi0 = reshape(eye(6), [36,1]);
    state0 = [V_soln(1:6); phi0];
    [tout, state_out] = ode113(@(t, state)CR3BP_full(state, mu), [0, V_soln(end)], state0, options);
    statef = state_out(end, :)';
    state_dot = CR3BP_full(statef, mu);
    DF_d = DF_mat(statef, state_dot);
    prev_n_hat = null(DF_d);

    family_member = 1;
    max_family_members = 700;

    repeat = true;

    while ((~p1_or_p2_bound) && (family_member <= max_family_members))
        print_out = sprintf('Family member - %d', family_member);
        disp(print_out)
        Vd = V_family(:,family_member);
        state0 = [Vd(1:6); phi0];
        [tout, state_out] = ode113(@(t, state)CR3BP_full(state, mu), [0, Vd(end)], state0, options);
        statef = state_out(end, :)';
        state_dot = CR3BP_full(statef, mu);
        DF_d = DF_mat(statef, state_dot);
        n_hat = null(DF_d);

        if dot(n_hat, prev_n_hat) < 0
            n_hat = -1 * n_hat;
        end

        n_hat = n_hat/norm(n_hat);

        V_star = V_family(:,family_member);
        V = V_star + delta_s*n_hat;
        H_norm = norm(H_psal(V_family(:,family_member), V_star, delta_s, n_hat, statef));
        counter = 1;
        
        repeat = true;
        runs = 1;
        while (repeat && runs < 10)
            try
                V = continuation(H_norm, TOL, counter, V, phi0, options, mu, n_hat, V_star, delta_s);
                repeat = false;
            catch
                delta_s = delta_s/2;
                V = V_star + delta_s*n_hat;
                repeat = true;
            end
            runs = runs + 1;
        end

        prev_n_hat = n_hat;

        V_family(:,family_member+1) = V(:,end);

        % Check if initial position is outside the P1, P2 bodies' normalized radius
        [tout, xout] = ode113(@(t, state)CR3BP(state, mu), [0, V_family(end,family_member+1)], V_family(1:6, family_member+1), options);
        for i = 1:length(tout)
            p2_minus_pos = p2_pos' - xout(i,1:3);
            if (norm(p2_minus_pos) < r_Moon_normalized)
                p1_or_p2_bound = true;
                disp("Hitting the Moon's surface.")
            end
            p1_minus_pos = p1_pos' - xout(i,1:3);
            if (norm(p1_minus_pos) < r_E_normalized)
                p1_or_p2_bound = true;
                disp("Hitting the Earth's surface.")
            end
        end

        family_member = family_member + 1;
    end
end

function V = continuation(H_norm, TOL, counter, V, phi0, options, mu, n_hat, V_star, delta_s)
    while ((H_norm > 1e-10) && (counter < 20))
        V_i = V(:,counter);
        state0 = [V_i(1:6); phi0];

        [tout, state_out] = ode113(@(t, state)CR3BP_full(state, mu), [0, V_i(end)], state0, options);
        
        statef = state_out(end, :)';
        state_dot = CR3BP_full(statef, mu);
        H_V = H_psal(V_i, V_star, delta_s, n_hat, statef);
        H_norm = norm(H_V);
        DH_V = DH_mat_psal(statef, state_dot, n_hat);
        V_ip1 = V_i - inv(DH_V) * H_V;
        V(:,counter+1) = V_ip1;
        counter = counter + 1;
    end
end

