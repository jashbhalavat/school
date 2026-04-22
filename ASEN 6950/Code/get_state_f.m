function statef_V = get_state_f(V, V_config, options, system_params)
    % Get statef for a given V

    mu = system_params(1);
    f = system_params(7);
    mdot = system_params(8);

    % V_config guides the number of arcs and whether each arc is thrusting
    % or not
    num_arcs = length(V_config);

    % Initialize statef_V
    % statef_V = [];
    statef_V = cell(num_arcs,1);
    
    % Start index at 1 to keep track of starting index for each arc
    % start_index = 1;
    for i = 1:length(V_config)
        if V_config(i) == 0
            % Natural arcs have 8 free variables - x, m, t
            % V_i = V(start_index:start_index+7);
            % start_index = start_index + 8;
            V_i = V{i};
            x_i_0 = V_i(1:6);
            m_i_0 = V_i(7);
            Dti = V_i(8);

            % Now propagate ICs
            [~, xi] = ode113(@(t,state)CR3BP(state, mu), [0, Dti], x_i_0, options);
            x_i_f = xi(end,1:6)';

            % Now populate statef
            statef_V{i} = [x_i_f; m_i_0; Dti];
        
        elseif V_config(i) == 1
            % Thrust arcs have 11 free variables - x, m, uhat, t
            % V_i = V(start_index:start_index+10);
            % start_index = start_index + 11;
            V_i = V{i};
            x_i_0 = V_i(1:6);
            m_i_0 = V_i(7);
            uhat_i_0 = V_i(8:10);
            Dti = V_i(11);

            % Dynamics function
            fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, uhat_i_0, f, mdot);

            % Propagate and get final value
            [~, xi] = ode113(fun, [0, Dti], [x_i_0; m_i_0], options);
            x_i_f = xi(end,1:6)';
            m_i_f = xi(end,7);

            % Now populate statef
            statef_V{i} = [x_i_f; m_i_f; uhat_i_0; Dti];
        end
    end


end