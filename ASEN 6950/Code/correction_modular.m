function V_soln = correction_modular(V0, V_config, system_params, x_init_des, x_final_des)
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
    f = system_params(7);
    mdot = system_params(8);

    % Set tolerance for numerical integrator and constraint vector
    TOL = 1e-12;

    % Set options for ode113
    options = odeset('RelTol', TOL, 'AbsTol', TOL);

    % Matrix of all free variable vectors
    % V(:,1) = V0;

    % While loop params
    counter = 1;
    counter_max = 50;

    % Get final state
    statef_V0 = get_state_f(V0, V_config, options, system_params);

    V = cell(counter_max,1);
    V{1} = V0;
    
    F_norm(1) = norm(cell2vec(F_modular(V0, statef_V0, V_config, x_init_des, x_final_des)));

    TOL = 1e-10;

    % While loop to reduce F_norm
    while ((F_norm(counter) > TOL) && (counter < counter_max))
        disp("Counter is " + counter);
        
        statef_V = get_state_f(V{counter}, V_config, options, system_params);

        F_i = F_modular(V{counter}, statef_V, V_config, x_init_des, x_final_des);
        DF_i = cell2mat(DF_mat_modular(V{counter}, V_config, F_i, options, system_params));

        V_i_vec = cell2vec(V{counter});
        F_i_vec = cell2vec(F_i);
        V_out = V_i_vec - DF_i' * inv(DF_i * DF_i') * F_i_vec;
        V{counter+1} = vec2cell(V_out, V_config);
    
        % Calculate F_norm and update counter
        F_norm(counter+1) = norm(F_i_vec);
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
