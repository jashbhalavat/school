function V_soln = correction(V0, system_params, T, Isp, uncorrected_init_thrust, uncorrected_final_thrust, t_star_em, l_star_em, init_mass, uncorrected_final_mass, desired_stage_1, desired_stage_4)
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
    final_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uncorrected_final_thrust, t_star_em, l_star_em, init_mass);

    % Propagate V0 non-linear CR3BP EOMs
    [~, xl1] = ode113(@(t,state)CR3BP(state, mu), [0, V0(7)], V0(1:6), options);
    [~, x1] = ode113(init_fun, [0, V0(14)], [V0(8:13); init_mass], options);
    [~, x2] = ode113(final_fun, [0, V0(21)], [V0(15:20); uncorrected_final_mass], options);
    [~, xl2] = ode113(@(t,state)CR3BP(state, mu), [0, V0(28)], V0(22:27), options);
    
    % Final final variables using V0
    statef_xl1 = xl1(end,1:6)';
    statef_x1 = x1(end,1:6)';
    statef_x2 = x2(end,1:6)';
    statef_xl2 = xl2(end,1:6)';
    statef_V0 = [statef_xl1; V0(7); statef_x1; V0(14); statef_x2; V0(21); statef_xl2; V0(28)];

    % Period is a free variable
    % T = V0(end);

    % Initialize constraint vector norm
    F_norm(1) = norm(F(V0, statef_V0, desired_stage_1, desired_stage_4));
    
    % Matrix of all free variable vectors
    V(:,1) = V0;

    % While loop params
    counter = 1;
    counter_max = 10;

    phi0 = reshape(eye(6), [36, 1]); % Initial phi is identity

    TOL = 1e-10;

    % While loop to reduce F_norm
    while ((F_norm(counter) > TOL) && (counter < counter_max))

        % Propagate full state and STM
        [~, xl1_out] = ode113(@(t,state)CR3BP(state, mu), [0, V(7,counter)], V(1:6,counter), options);
        [~, x1_out] = ode113(init_fun, [0 V(14,counter)], [V(8:13,counter); init_mass], options);
        [~, x2_out] = ode113(final_fun, [0 V(21,counter)], [V(15:20,counter); uncorrected_final_mass], options);
        [~, xl2_out] = ode113(@(t,state)CR3BP(state, mu), [0, V(28,counter)], V(22:27,counter), options);

        xl1_f = xl1_out(end, 1:6)';
        x1_f = x1_out(end, 1:6)';
        x2_f = x2_out(end, 1:6)';
        xl2_f = xl2_out(end, 1:6)';
        statef = [xl1_f; V(7,counter);
            x1_f; V(14,counter);
            x2_f; V(21,counter);
            xl2_f; V(28,counter);];

        F_i = F(V(:,counter), statef, desired_stage_1, desired_stage_4);
        DF_i = DF_mat(V(:,counter), options, mu, T, Isp, uncorrected_init_thrust, uncorrected_final_thrust, t_star_em, l_star_em, init_mass, uncorrected_final_mass, desired_stage_1, desired_stage_4);

        % Find V_i+1
        V(:,counter+1) = V(:,counter) - DF_i' * inv(DF_i * DF_i') * F_i;
    
        % Calculate F_norm and update counter
        F_norm(counter+1) = norm(F_i);
        counter = counter + 1;
    end

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

    V_soln = V(:,end);

end
