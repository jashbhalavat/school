function plot_full_comparison_modular(V0, V_soln, V_config, l1_pos, l2_pos, system_params)
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

    xout = cell(length(V_config),2);
    tout = cell(length(V_config),2);
    m = cell(length(V_config),2);
    uhat = cell(length(V_config),2);

    for i = 1:length(V_config)
        % Go through all the free variables and plot the initial guess
        V_post = V_soln{i};
        V_pre = V0{i};
        if V_config(i) == 0
            % Natural arc
            [tout_i_pre, xout{i,1}] = ode113(@(t,state)CR3BP(state, mu), [0, V_pre(8)], V_pre(1:6), options_no_events);
            [tout_i_post, xout{i,2}] = ode113(@(t,state)CR3BP(state, mu), [0, V_post(8)], V_post(1:6), options_no_events);
            uhat{i,1} = zeros([length(tout_i_pre), 3]);
            uhat{i,2} = zeros([length(tout_i_post), 3]);
            if i == 1
                m{i,1} = ones(size(tout_i_pre));
                tout{i,1} = tout_i_pre;
                m{i,2} = ones(size(tout_i_post));
                tout{i,2} = tout_i_post;
            else
                last_mass_pre = m{i-1,1}(end);
                m{i,1} = last_mass_pre*ones(size(tout_i_pre));
                last_time_pre = tout{i-1,1}(end);
                tout{i,1} = tout_i_pre + last_time_pre;
                last_mass_post = m{i-1,2}(end);
                m{i,2} = last_mass_post*ones(size(tout_i_post));
                last_time_post = tout{i-1,2}(end);
                tout{i,2} = tout_i_post + last_time_post;
            end
        elseif V_config(i) == 1
            % Thrusting arc
            fun_pre = @(t,state)CR3BP_with_non_dim_mass(state, mu, V_pre(8:10), f, mdot);
            [tout_i_pre, xout_i_pre] = ode113(fun_pre, [0, V_pre(11)], V_pre(1:7), options_no_events);
            m{i,1} = xout_i_pre(:,7);
            xout{i,1} = xout_i_pre(:,1:6);
            uhat{i,1} = V_pre(8:10)'.*ones([length(tout_i_pre), 3]);

            fun_post = @(t,state)CR3BP_with_non_dim_mass(state, mu, V_post(8:10), f, mdot);
            [tout_i_post, xout_i_post] = ode113(fun_post, [0, V_post(11)], V_post(1:7), options_no_events);
            m{i,2} = xout_i_post(:,7);
            xout{i,2} = xout_i_post(:,1:6);
            uhat{i,2} = V_post(8:10)'.*ones([length(tout_i_post), 3]);
            if i == 1
                tout{i,1} = tout_i_pre;
                tout{i,2} = tout_i_post;
            else
                last_time_pre = tout{i-1,1}(end);
                tout{i,1} = tout_i_pre + last_time_pre;
                last_time_post = tout{i-1,2}(end);
                tout{i,2} = tout_i_post + last_time_post;
            end
        end
    end

    xout_pre_vec = cell2vec(xout(:,1));
    tout_pre_vec = cell2vec(tout(:,1));
    mass_pre_vec = cell2vec(m(:,1));
    uhat_pre_vec = cell2vec(uhat(:,1));

    xout_post_vec = cell2vec(xout(:,2));
    tout_post_vec = cell2vec(tout(:,2));
    mass_post_vec = cell2vec(m(:,2));
    uhat_post_vec = cell2vec(uhat(:,2));

    % First figure is the trajectory
    figure()
    scatter(l1_pos(1), l1_pos(2), 'filled', 'cyan')
    hold on
    scatter(l2_pos(1), l2_pos(2), 'filled', 'black')
    % plot(init_orbit(:,1), init_orbit(:,2), 'blue', 'LineWidth', 2)
    % plot(final_orbit(:,1), final_orbit(:,2), 'red', 'LineWidth', 2)
    grid on
    axis equal
    title("Transfer")
    xlabel('$$\hat{x}$$','Interpreter','Latex', 'FontSize',18)
    ylabel('$$\hat{y}$$','Interpreter','Latex', 'FontSize',18)
    plot(xout_pre_vec(:,1), xout_pre_vec(:,2),'LineStyle', '--', 'color','red','LineWidth', 2)
    plot(xout_post_vec(:,1), xout_post_vec(:,2), 'blue','LineWidth', 1)
    legend("L1", "L2", "Uncorrected Transfer", "Corrected Transfer")

    % Mass
    figure()
    plot(tout_pre_vec, mass_pre_vec,'LineStyle', '--', 'color','red','LineWidth', 2)
    hold on
    plot(tout_post_vec, mass_post_vec, 'blue', 'LineWidth',1)
    grid on
    legend("Uncorrected Transfer", "Corrected Transfer")
    xlabel("Non-dimenstional time")
    ylabel("Normalized mass")
    title("Mass vs Time")

    % Jacobian
    jacobian_pre = jacobiConstantCR3BP(xout_pre_vec, mu);
    jacobian_post = jacobiConstantCR3BP(xout_post_vec, mu);
    figure()
    plot(tout_pre_vec, jacobian_pre,'LineStyle', '--', 'color','red','LineWidth', 2)
    hold on
    plot(tout_post_vec, jacobian_post, 'blue', 'LineWidth',1)
    grid on
    xlabel("Non-dimenstional time")
    legend("Uncorrected Transfer", "Corrected Transfer")
    ylabel("Jacobian")
    title("Jacobian vs Time")

    % Control Vector
    figure()
    subplot(3,1,1)
    plot(tout_pre_vec, uhat_pre_vec(:,1),'LineStyle', '--', 'color','red','LineWidth', 2)
    hold on
    plot(tout_post_vec, uhat_post_vec(:,1),'LineStyle', '--', 'color','blue','LineWidth', 1)
    grid on
    xlabel("Non-dimenstional time")
    ylabel("uHat_x")
    legend("Uncorrected Transfer", "Corrected Transfer")
    sgtitle("Thrust Vector vs Time")

    subplot(3,1,2)
    plot(tout_pre_vec, uhat_pre_vec(:,2),'LineStyle', '--', 'color','red','LineWidth', 2)
    hold on
    plot(tout_post_vec, uhat_post_vec(:,2),'LineStyle', '--', 'color','blue','LineWidth', 1)
    grid on
    xlabel("Non-dimenstional time")
    ylabel("uHat_y")

    subplot(3,1,3)
    plot(tout_pre_vec, uhat_pre_vec(:,3),'LineStyle', '--', 'color','red','LineWidth', 2)
    hold on
    plot(tout_post_vec, uhat_post_vec(:,3),'LineStyle', '--', 'color','blue','LineWidth', 1)
    grid on
    xlabel("Non-dimenstional time")
    ylabel("uHat_z")




end