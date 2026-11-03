function plot_full_modular(V, V_config, l1_pos, l2_pos, system_params)
    % Function to plot initial guess and corrected trajectory along with
    % mass, thrust vector, and jacobi constant

    % system_params = [mu, t_star_em, l_star_em, T, Isp, init_mass, f, mdot];
    mu = system_params(1);
    t_star_em = system_params(2);
    l_star_em = system_params(3);
    T = system_params(4);
    Isp = system_params(5);
    init_mass = system_params(6);
    f = system_params(7);
    mdot = system_params(8);

    options_no_events = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    xout = cell(length(V_config),1);
    tout = cell(length(V_config),1);
    m = cell(length(V_config),1);

    for i = 1:length(V_config)
        % Go through all the free variables and plot the initial guess
        V_curr = V{i};
        if V_config(i) == 0
            % Natural arc
            [tout_i, xout{i}] = ode113(@(t,state)CR3BP(state, mu), [0, V_curr(8)], V_curr(1:6), options_no_events);
            if i == 1
                m{i} = ones(size(tout_i));
                tout{i} = tout_i;
            else
                last_mass = m{i-1}(end);
                m{i} = last_mass*ones(size(tout_i));
                last_time = tout{i-1}(end);
                tout{i} = tout_i + last_time;
            end
            % plot(xout(:,1), xout(:,2),'blue', 'LineWidth', 2)
        elseif V_config(i) == 1
            % Thrusting arc
            fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, V_curr(8:10), f, mdot);
            [tout_i, xout_i] = ode113(fun, [0, V_curr(11)], V_curr(1:7), options_no_events);
            m{i} = xout_i(:,7);
            xout{i} = xout_i(:,1:6);
            if i == 1
                tout{i} = tout_i;
            else
                last_time = tout{i-1}(end);
                tout{i} = tout_i + last_time;
            end
            % plot(xout(:,1), xout(:,2),'red' ,'LineWidth', 2)
        end
    end

    xout_vec = cell2vec(xout);
    tout_vec = cell2vec(tout);
    mass_vec = cell2vec(m);

    % First figure is the trajectory
    figure()
    scatter(l1_pos(1), l1_pos(2), 'filled', 'red')
    hold on
    scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
    % plot(init_orbit(:,1), init_orbit(:,2), 'blue', 'LineWidth', 2)
    % plot(final_orbit(:,1), final_orbit(:,2), 'red', 'LineWidth', 2)
    grid on
    axis equal
    title("Transfer")
    xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
    plot(xout_vec(:,1), xout_vec(:,2), 'blue', 'LineWidth', 2)
    legend("L1", "L2", "Transfer")

    % Mass
    figure()
    plot(tout_vec, mass_vec, 'LineWidth',2)
    grid on
    xlabel("Non-dimenstional time")
    ylabel("Normalized mass")

    % Jacobian
    jacobian = jacobiConstantCR3BP(xout_vec, mu);
    figure()
    plot(tout_vec, jacobian, 'LineWidth',2)
    grid on
    xlabel("Non-dimenstional time")
    ylabel("Jacobian")





end