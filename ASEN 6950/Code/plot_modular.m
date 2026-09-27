function plot_modular(V, V_config, l1_pos, l2_pos, system_params)
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

    % First figure is the trajectory
    figure()
    scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
    hold on
    scatter(l1_pos(1), l1_pos(2), 'filled', 'black')
    % plot(init_orbit(:,1), init_orbit(:,2), 'blue', 'LineWidth', 2)
    % plot(final_orbit(:,1), final_orbit(:,2), 'red', 'LineWidth', 2)
    grid on
    axis equal
    title("Transfer")
    xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)

    for i = 1:length(V_config)
        % Go through all the free variables and plot the initial guess
        V_curr = V{i};
        if V_config(i) == 0
            % Natural arc
            [tout, xout] = ode113(@(t,state)CR3BP(state, mu), [0, V_curr(8)], V_curr(1:6), options_no_events);
            plot(xout(:,1), xout(:,2), 'LineWidth', 2)
        elseif V_config(i) == 1
            % Thrusting arc
            fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, V_curr(8:10), f, mdot);
            [tout, xout] = ode113(fun, [0, V_curr(11)], V_curr(1:7), options_no_events);
            plot(xout(:,1), xout(:,2), 'LineWidth', 2)
        end
    end

end