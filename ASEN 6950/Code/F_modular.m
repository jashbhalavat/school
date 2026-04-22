function F = F_modular(state0, statef, V_config, x_init_des, x_final_des)
    % Modified Constraint Vector

    num_arcs = length(V_config);

    F = cell(num_arcs+1,1);

    F{1} = [state0{1}(1:6) - x_init_des; % Anchor initial state
            state0{1}(7) - 1]; % Anchor initial mass

    for i = 1:num_arcs
        if i == num_arcs
            % If final arc, anchor to final point
            
            % Get current arc free variables
            % V_i_0 = state0(start_index:end);
            % V_i_f = statef(start_index:end);
            V_i_0 = state0{i};
            V_i_f = statef{i};

            % For natural arcs, we're only constraining the state
            F_i = V_i_f(1:6) - x_final_des;

            if V_config(i) == 1
                % If this is a thrusting arc, also constrain uhat magnitude
                % NOTE - not making sure final mass is >1 as of
                % 20260415
                F_i = [F_i; norm(V_i_0(8:10))^2 - 1];
            end
            F{i+1} = F_i;
        else
            % If not a final arc, constrain between arcs

            if V_config(i) == 0
                % For natural arcs only constraining x and m

                % V_i_f = statef(start_index:start_index+7);
                % start_index = start_index + 8;
                V_i_f = statef{i};
                x_i_f = V_i_f(1:6);
                m_i_f = V_i_f(7);

                V_ip1_0 = state0{i+1};
                x_ip1_0 = V_ip1_0(1:6);
                m_ip1_0 = V_ip1_0(7);

                F_i = [x_i_f - x_ip1_0;
                        m_i_f - m_ip1_0];
                F{i+1} = F_i;

            elseif V_config(i) == 1
                % For thrust arcs constraining x, m, and u_mag

                % V_i_f = statef(start_index:start_index+10);
                % start_index = start_index + 11;
                V_i_f = statef{i};
                x_i_f = V_i_f(1:6);
                m_i_f = V_i_f(7);
                uhat_i = V_i_f(8:10);
                
                V_ip1_0 = state0{i+1};
                x_ip1_0 = V_ip1_0(1:6);
                m_ip1_0 = V_ip1_0(7);

                F_i = [x_i_f - x_ip1_0;
                        m_i_f - m_ip1_0;
                        norm(uhat_i)^2 - 1];
                F{i+1} = F_i;
            end
        end
    end

end