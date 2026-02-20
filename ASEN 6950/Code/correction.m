function V_soln = correction(V0, system_params, T, Isp, uncorrected_init_thrust, uncorrected_final_thrust, t_star_em, l_star_em, init_mass, uncorrected_final_mass)
    % Script to compute a general three-dimensional periodic orbit via multiple shooting
    % Inputs
    % V0 - initial guess for a free variable vector - [x_1_0, delta_t1,
    % x_2_0, delta_t_2]

    % system_params - system parameters
    % 
    % Output
    % V_soln - free variable vector corresponding to a solution

    % Get mass ratio of system
    mu = system_params(1);

    % Set tolerance for numerical integrator and constraint vector
    TOL = 1e-12;

    % Set options for ode113
    options = odeset('RelTol', TOL, 'AbsTol', TOL);

    init_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_init_thrust, t_star_em, l_star_em, init_mass);
    final_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, uncorrected_final_mass);

    % Propagate V0 non-linear CR3BP EOMs
    [~, x1] = ode113(init_fun, [0, V0(8)], V0(1:7), options);
    [~, x2] = ode113(final_fun, [0, V0(16)], V0(9:15), options);
    
    % Final final variables using V0
    statef_x1 = x1(end,:)';
    statef_x2 = x2(end,:)';
    statef_V0 = [statef_x1; V0(8); statef_x2; V0(16)];

    % Period is a free variable
    % T = V0(end);

    % Initialize constraint vector norm
    F_norm(1) = norm(F(V0, statef_V0));
    
    % Matrix of all free variable vectors
    V(:,1) = V0;

    % While loop params
    counter = 1;
    counter_max = 50;

    phi0 = reshape(eye(6), [36, 1]); % Initial phi is identity

    TOL = 1e-10;

    % While loop to reduce F_norm
    while ((F_norm(counter) > TOL) && (counter < counter_max))

        init_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_init_thrust, t_star_em, l_star_em, V(7,counter));
        final_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, V(15,counter));
        
        % Propagate full state and STM
        [~, x1_out] = ode113(init_fun, [0 V(8,counter)], V(1:7,counter), options);
        [~, x2_out] = ode113(final_fun, [0 V(16,counter)], V(9:15,counter), options);

        x1_f = x1_out(end, :)';
        x2_f = x2_out(end, :)';
        statef = [x1_f; V(8,counter); x2_f; V(16,counter)];

        F_i = F(V(:,counter), statef);
        DF_i = DF_mat(V(:,counter), options, mu);

        % Find V_i+1
        V(:,counter+1) = V(:,counter) - DF_i' * inv(DF_i * DF_i') * F_i;
    
        % Calculate F_norm and update counter
        F_norm(counter+1) = norm(F_i);
        counter = counter + 1;
    end

    if plot_input
        figure()
        plot([1:counter], F_norm, '-o', 'LineWidth', 2)
        yscale log
        grid on
        xlabel("Iterations")
        ylabel("F Norm")
        title("Constraint Vector Norm for each Iteration")
        hold on
        tol_yline = ones([counter,1])*TOL;

        plot([1:counter], tol_yline, 'red', 'LineWidth', 2)
        hold off
        legend("Norm", "Threshold")
    end

    V_soln = V(:,end);

end
