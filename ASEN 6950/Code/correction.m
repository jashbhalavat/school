% function V_soln = correction(V0, system_params, T, Isp, uncorrected_init_thrust, uncorrected_final_thrust, t_star_em, l_star, init_mass, uncorrected_final_mass, desired_stage_1, desired_stage_4)
% function V_soln = correction(V0, system_params, x_1_des, x_4_des)
function V_soln = correction(V0, system_params, x_1_des, x_2_des, x_3_des)
    % Script to compute a general three-dimensional periodic orbit via multiple shooting
    % Inputs
    % V0 - initial guess for a free variable vector
    %       [x_1,0 (1:6); Dt1 (7); x_2,0 (8:13); m_2,0 (14); Dt2 (15); uhat_2 (16:18); x_3,0 (19:24); m_3,0 (25); Dt3 (26);
    %       uhat_3 (27:29); x_4,0 (30:35); Dt4 (36)]

    % system_params - mu, t_star, l_star, T, Isp, init_mass
    % 
    % Output
    % V_soln - free variable vector corresponding to a solution

    % Get mass ratio of system
    mu = system_params(1);
    t_star = system_params(2);
    l_star = system_params(3);
    T = system_params(4);
    Isp = system_params(5);
    M_sc_0 = system_params(6);

    % Set tolerance for numerical integrator and constraint vector
    TOL = 1e-12;

    % Set options for ode113
    options = odeset('RelTol', TOL, 'AbsTol', TOL);

    % Matrix of all free variable vectors
    V(:,1) = V0;

    x_1_0 = V0(1:6);
    Dt1 = V0(7);
    x_2_0 = V0(8:13);
    m_2_0 = V0(14);
    Dt2 = V0(15);
    uhat_2 = V0(16:18);
    x_3_0 = V0(19:24);
    m_3_0 = V0(25);
    Dt3 = V0(26);
    uhat_3 = V0(27:29);
    % x_4_0 = V0(30:35);
    % Dt4 = V0(36);

    % Define functions
    init_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uhat_2, t_star, l_star, m_2_0);
    final_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uhat_3, t_star, l_star, m_3_0);

    % Propagate V0 non-linear CR3BP EOMs
    [~, x1] = ode113(@(t,state)CR3BP(state, mu), [0, Dt1], x_1_0, options);
    [~, x2] = ode113(init_fun, [0, Dt2], [x_2_0; m_2_0], options);
    [~, x3] = ode113(final_fun, [0, Dt3], [x_3_0; m_3_0], options);
    % [~, x4] = ode113(@(t,state)CR3BP(state, mu), [0, Dt4], x_4_0, options);

    % Final states
    x_1_f = x1(end,1:6)';
    x_2_f = x2(end,1:6)';
    m_2_f = x2(end,7);
    x_3_f = x3(end,1:6)';
    m_3_f = x3(end,7);
    % x_4_f = x4(end,1:6)';
    
    % statef_V0 = [x_1_f; Dt1; x_2_f; m_2_f; Dt2; uhat_2; x_3_f; m_3_f; Dt3; uhat_3; x_4_f; Dt4];
    statef_V0 = [x_1_f; Dt1; x_2_f; m_2_f; Dt2; uhat_2; x_3_f; m_3_f; Dt3; uhat_3];

    % While loop params
    counter = 1;
    counter_max = 50;
    
    % F_norm(1) = norm(F(V0, statef_V0, x_1_des, x_4_des));
    F_norm(1) = norm(F(V0, statef_V0, x_1_des, x_2_des, x_3_des, M_sc_0));

    phi0 = reshape(eye(6), [36, 1]); % Initial phi is identity

    TOL = 1e-10;

    % While loop to reduce F_norm
    while ((F_norm(counter) > TOL) && (counter < counter_max))
        disp("Counter is " + counter);
        x_1_0 = V(1:6,counter);
        Dt1 = V(7,counter);
        x_2_0 = V(8:13,counter);
        m_2_0 = V(14,counter);
        Dt2 = V(15,counter);
        uhat_2 = V(16:18,counter);
        x_3_0 = V(19:24,counter);
        m_3_0 = V(25,counter);
        Dt3 = V(26,counter);
        uhat_3 = V(27:29,counter);
        % x_4_0 = V(30:35,counter);
        % Dt4 = V(36,counter);

        % Define functions
        init_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uhat_2, t_star, l_star, m_2_0);
        final_fun = @(t,state)CR3BP_with_non_dim_mass(state, mu, T, Isp, uhat_3, t_star, l_star, m_3_0);

        % Propagate V0 non-linear CR3BP EOMs
        [~, x1] = ode113(@(t,state)CR3BP(state, mu), [0, Dt1], x_1_0, options);
        [~, x2] = ode113(init_fun, [0, Dt2], [x_2_0; m_2_0], options);
        [~, x3] = ode113(final_fun, [0, Dt3], [x_3_0; m_3_0], options);
        % [~, x4] = ode113(@(t,state)CR3BP(state, mu), [0, Dt4], x_4_0, options);

        % Final states
        x_1_f = x1(end,1:6)';
        x_2_f = x2(end,1:6)';
        m_2_f = x2(end,7);
        x_3_f = x3(end,1:6)';
        m_3_f = x3(end,7);
        % x_4_f = x4(end,1:6)';

        % statef_V = [x_1_f; Dt1; x_2_f; m_2_f; Dt2; uhat_2; x_3_f; m_3_f; Dt3; uhat_3; x_4_f; Dt4];
        statef_V = [x_1_f; Dt1; x_2_f; m_2_f; Dt2; uhat_2; x_3_f; m_3_f; Dt3; uhat_3];

        % F_i = F(V(:,counter), statef_V, x_1_des, x_4_des);
        % DF_i = DF_mat(V(:,counter), statef_V, options, system_params, x_1_des, x_4_des);

        F_i = F(V(:,counter), statef_V, x_1_des, x_2_des, x_3_des, M_sc_0);
        DF_i = DF_mat(V(:,counter), statef_V, options, system_params, x_1_des, x_2_des, x_3_des);

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
    tol_yline = ones([counter+5,1])*TOL;

    plot([1:counter+5], tol_yline, 'red', 'LineWidth', 2)
    hold off
    legend("Norm", "Threshold")

    V_soln = V(:,end);

end
