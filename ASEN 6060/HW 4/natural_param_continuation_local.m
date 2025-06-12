function V_family = natural_param_continuation_local(V_soln, mu)
    % V_soln - First corrected free variable solution
    % mu - system parameter

    % out - Family of free variables

    TOL = 1e-12;

    V_family = V_soln;
    delta = -1e-3;

    family_members = 50;

    % Set options for ode113
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    [tout_corrected, xout_corrected] = ode113(@(t, state)CR3BP(state, mu), [0, V_soln(end)], V_soln(1:6), options);

    for i = 1:family_members
        print_out = sprintf('Family member - %d', i);
        disp(print_out)
        V = V_family(:,i);
        x_star = V(1);
        x_star_plus_delta = x_star + delta;
        H_norm = norm(H(V_soln(1:6), xout_corrected(end,:), x_star_plus_delta));
        counter = 1;
        
        while ((H_norm > TOL) && (counter < 20))
            V_i = V(:,counter);
            phi_0 = reshape(eye(6), [36,1]);
            state0 = [V_i(1:6); phi_0];

            [tout, state_out] = ode113(@(t, state)CR3BP_full(state, mu), [0, V_i(end)], state0, options);
            
            statef = state_out(end, :)';
            state_dot = CR3BP_full(statef, mu);
            H_V = H(V_i(1:6), statef, x_star_plus_delta);
            H_norm = norm(H_V);
            DH_V = DH_mat(statef, state_dot);
            V_ip1 = V_i - inv(DH_V) * H_V;
            V(:,counter+1) = V_ip1;
            counter = counter + 1;
        end

        V_family(:,i+1) = V(:,end);
    end
end
    