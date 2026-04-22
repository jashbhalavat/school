function DF = DF_mat_modular(V, V_config, F_i, options, system_params)
    % Modified constraint DF matrix

    % Get mass ratio of system
    mu = system_params(1);
    t_star = system_params(2);
    l_star = system_params(3);
    T = system_params(4);
    Isp = system_params(5);
    M_sc_0 = system_params(6);
    f = system_params(7);
    mdot = system_params(8);

    num_arcs = length(V_config);

    phi0_36 = reshape(eye(6), [36, 1]); % Initial phi is identity
    phi0_100 = reshape(eye(10), [100, 1]); % Initial phi is identity

    xdot_f = cell(num_arcs,1);
    phi_mat = cell(num_arcs,1);

    for i = 1:num_arcs
        if V_config(i) == 0
            % Natural arcs have 8 free variables - x, m, t
            % V_i = V(start_index:start_index+7);
            % start_index = start_index + 8;
            V_i = V{i};
            x_i_0 = V_i(1:6);
            m_i_0 = V_i(7);
            Dti = V_i(8);

            xi_0 = [x_i_0; phi0_36];
            [~, xi_out] = ode113(@(t,state)CR3BP_full(state, mu), [0, Dti], xi_0, options);
            xi_f = xi_out(end, :)';
            xdot_i_f = CR3BP(xi_f, mu);
            xdot_f{i} = xdot_i_f;

            phi_row_i = xi_f(7:end);
            phi_mat{i} = reshape(phi_row_i, [6,6])';

        elseif V_config(i) == 1
            % Thrust arcs have 11 free variables - x, m, uhat, t
            % V_i = V(start_index:start_index+10);
            % start_index = start_index + 11;
            V_i = V{i};
            x_i_0 = V_i(1:6);
            m_i_0 = V_i(7);
            uhat_i_0 = V_i(8:10);
            Dti = V_i(11);

            xi_0 = [x_i_0; m_i_0; uhat_i_0; phi0_100];
            fun = @(t,state)CR3BP_full_mass(state, mu, T, Isp, t_star, l_star, M_sc_0);
            [~, xi_out] = ode113(fun, [0 Dti], xi_0, options);
            xi_f = xi_out(end, :)';
            xdot_i_f = CR3BP_with_non_dim_mass(xi_f, mu, uhat_i_0, f, mdot);
            xdot_f{i} = xdot_i_f;

            phi_row_i = xi_f(11:end);
            phi_mat{i} = reshape(phi_row_i, [10,10])';
        end
    end


    % Instantiate DF matrix
    % DF = zeros(sum(cellfun(@numel,F_i)), sum(cellfun(@numel,V)));
    DF = cell(length(F_i), length(V));

    for i = 1:length(F_i)
        % F_i length will always be one more than V (which should have
        % length of num_arcs)
        for j = 1:length(V)
            DF{i, j} = zeros(length(F_i{i}), length(V{j}));
        end
    end

    % The first 7 diags are always 1 because of anchoring to initial state
    % and mass
    DF{1,1}(1:7,1:7) = eye(7);

    % Now, go through all the individual blocks and fill out based on
    % whether the arc is thrusting or natural
    for i = 2:length(F_i)
        % start at 2 since the first block of F_i was identity
        if i == num_arcs+1
            % Last constraint vector
            if V_config(i-1) == 0
                % Natural arc
                DF{i,i-1} = [phi_mat{i-1}, zeros(6,1), xdot_f{i-1}];
            elseif V_config(i-1) == 1
                % Thrusting arc
                uhat_i = V{i-1}(8:10);
                DF{i,i-1} = [phi_mat{i-1}(1:6,:), xdot_f{i-1}(1:6);
                            zeros(1,7), 2*uhat_i', 0];
            end
        else
            % Not last constraint vector
            if V_config(i-1) == 0
                % Natural arc
                DF{i,i-1} = [phi_mat{i-1}, zeros(6,1), xdot_f{i-1};
                        zeros(1,6), 1, 0];
                DF{i,i}(1:7,1:7) = -eye(7);
            elseif V_config(i-1) == 1
                % Thrusting arc
                uhat_i = V{i-1}(8:10);
                DF{i,i-1} = [phi_mat{i-1}(1:7,:), xdot_f{i-1}(1:7);
                            zeros(1,7), 2*uhat_i', 0];
                DF{i,i}(1:7,1:7) = -eye(7);
            end
        end
    end

end