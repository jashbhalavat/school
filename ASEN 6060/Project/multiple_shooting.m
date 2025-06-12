function V_soln = multiple_shooting(V0, system_params, plot_input, r_des_4)
    % Script to compute a general three-dimensional periodic orbit via multiple shooting
    % Inputs
    % V0 - initial guess for a free variable vector
    % statef_V0 - final state when V0 is used as initial guess using CR3BRP
    % EOMs
    % system_params - system parameters
    % r_des_4 - Desired position vector after 4th arc
    % 
    % Output
    % V_soln - free variable vector corresponding to a solution

    % Get mass ratio of system
    mu = system_params(1);

    % Set tolerance for numerical integrator and constraint vector
    TOL = 1e-12;

    % Set options for ode113
    options = odeset('RelTol', TOL, 'AbsTol', TOL);

    % Propagate V0 non-linear CR3BP EOMs
    [tout, x1] = ode113(@(t, state)CR3BP(state, mu), [0 V0(7)], V0(1:6), options);
    [tout, x2] = ode113(@(t, state)CR3BP(state, mu), [0 V0(14)], V0(8:13), options);
    [tout, x3] = ode113(@(t, state)CR3BP(state, mu), [0 V0(21)], V0(15:20), options);
    [tout, x4] = ode113(@(t, state)CR3BP(state, mu), [0 V0(28)], V0(22:27), options);
    
    % Final final variables using V0
    % statef_V0 = xout(end,:);
    statef_x1 = x1(end,:)';
    statef_x2 = x2(end,:)';
    statef_x3 = x3(end,:)';
    statef_x4 = x4(end,:)';
    statef_V0 = [statef_x1; V0(7); statef_x2; V0(14); statef_x3; V0(21); statef_x4; V0(28)];

    % Period is a free variable
    % T = V0(end);

    % Initialize constraint vector norm
    F_norm(1) = norm(F(V0, statef_V0, r_des_4));
    
    % Matrix of all free variable vectors
    V(:,1) = V0;

    % While loop params
    counter = 1;
    counter_max = 50;

    phi0 = reshape(eye(6), [36, 1]); % Initial phi is identity

    TOL = 1e-10;

    % While loop to reduce F_norm
    while ((F_norm(counter) > TOL) && (counter < counter_max))
        % x1_0 = [V(1:6,counter); phi0];
        % x2_0 = [V(8:13,counter); phi0];
        % x3_0 = [V(15:20,counter); phi0];
        % x4_0 = [V(22:27,counter); phi0];
        
        % Propagate full state and STM
        [~, x1_out] = ode113(@(t, state)CR3BP(state, mu), [0 V(7,counter)], V(1:6,counter), options);
        [~, x2_out] = ode113(@(t,state)CR3BP(state, mu), [0 V(14,counter)], V(8:13,counter), options);
        [~, x3_out] = ode113(@(t,state)CR3BP(state, mu), [0 V(21,counter)], V(15:20,counter), options);
        [~, x4_out] = ode113(@(t,state)CR3BP(state, mu), [0 V(28,counter)], V(22:27,counter), options);

        x1_f = x1_out(end, :)';
        x2_f = x2_out(end, :)';
        x3_f = x3_out(end, :)';
        x4_f = x4_out(end, :)';
        statef = [x1_f; V(7,counter); x2_f; V(14,counter); x3_f; V(21,counter); x4_f; V(28,counter)];

        F_i = F(V(:,counter), statef, r_des_4);
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
