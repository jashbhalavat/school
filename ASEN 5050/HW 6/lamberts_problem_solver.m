function out = lamberts_problem_solver(mu, R1, V1, R2, V2, TOF, delta_a, greater_than_180)
    % Solve lambert's problem using fsolve
    % Inputs - gravitational parameter of central body,
    %          initial pos, vel, final pos, vel, Time of Flight,
    %          initial guess delta a, whether choice is greater/lower than
    %          180 deg
    % Output - Total delta V

    % Compute norms of initial and final positions
    r1 = norm(R1);
    r2 = norm(R2);

    % Step 1 - Transfer angle
    delta_theta_star = acos(dot(R1, R2) / (r1*r2));

    % If greater than 180, subtract original value by 360 deg to get a
    % delta theta star that's > 180 deg
    if greater_than_180
        delta_theta_star = 2*pi - delta_theta_star;
    end

    % Step 2 - calculate c and s
    c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(delta_theta_star));
    s = 0.5 * (r1 + r2 + c);

    TOF_P = 1/3 * sqrt(2/mu) * (s^(3/2) + (s - c)^(3/2));

    elliptical_orbit = false;
    if TOF > TOF_P
        elliptical_orbit = true;
    end

    if elliptical_orbit
        am = s/2;
        alpha_m = pi;
        nm = sqrt(mu/am^3);
        beta_m_0 = 2 * asin(sqrt((s-c)/s));
        beta_m = beta_m_0;
        if greater_than_180
            beta_m = -beta_m_0;
        end
        TOF_min = 1/nm * (alpha_m - beta_m - (sin(alpha_m) - sin(beta_m)));

        a_initial_guess = am + delta_a;

        fun = @(a)lamberts_eqn_elliptical(a, s, c, TOF, mu, TOF_min, greater_than_180);
        a = fsolve(fun, a_initial_guess);

        alpha_0_final = 2 * asin(sqrt(s/(2*a)));
        beta_0_final = 2 * asin(sqrt((s-c)/(2*a)));
        alpha_beta_final = lamberts_problem_alpha_beta_logic(TOF, TOF_min, greater_than_180, alpha_0_final, beta_0_final);
        alpha_final = alpha_beta_final(1);
        beta_final = alpha_beta_final(2);
        n_final = sqrt(mu/a^3);
        TOF_final = 1/n_final * (alpha_final - beta_final - (sin(alpha_final) - sin(beta_final)));
        if abs(TOF_final - TOF) > 1e-4
            disp("TOFs don't match\n")
        end
    else
        a_initial_guess = delta_a;
        fun = @(a)lamberts_eqn_hyperbolic(a, s, c, TOF, mu, greater_than_180);

        a = fsolve(fun, a_initial_guess);
        alpha_final = 2 * asinh(sqrt(s/(2*abs(a))));
        beta_final = 2 * asinh(sqrt((s-c)/(2*abs(a))));
        if greater_than_180
            TOF_final = sqrt(a^3/mu) * (sinh(alpha_final) - alpha_final + (sinh(beta_final) - beta_final));
        else
            TOF_final = sqrt(a^3/mu) * (sinh(alpha_final) - alpha_final - (sinh(beta_final) - beta_final));
        end

        if abs(TOF_final - TOF) > 1e-4
            disp("TOFs don't match\n")
        end
    end

    p = (4*a*(s-r1)*(s-r2))/(c^2) * (sin((alpha_final + beta_final)/2))^2;
    e = sqrt(1 - p/a);

    theta_star_1_p = abs(acos(1/e * (p/r1 - 1)));
    theta_star_2_p = abs(acos(1/e * (p/r2 - 1)));
    theta_star_1_n = -theta_star_1_p;
    theta_star_2_n = -theta_star_2_p;

    % Test matrix
    % First column = difference between 2 and 1
    % Second column = sign of theta_star_2
    % Third column = sign of theta_star_1
    test_ones = ones(6, 2);
    test = [zeros(6, 1), test_ones];
    test(1, 1) = theta_star_2_p - theta_star_1_p;
    test(2, 1) = theta_star_2_p - theta_star_1_n;
    test(2, 3) = -1;
    test(3, 1) = theta_star_2_n - theta_star_1_p;
    test(3, 2) = -1;
    test(4, 1) = theta_star_2_n - theta_star_1_n;
    test(4, 2) = -1;
    test(4, 3) = -1;
    test(5, 1) = 2*pi - theta_star_2_p - theta_star_1_p;
    test(5, 2) = -1;
    test(6, 1) = 2*pi - theta_star_2_p - theta_star_1_n;
    test(6, 2) = -1;
    test(6, 3) = -1;

    for i = 1:length(test)
        if abs(test(i) - delta_theta_star) < 1e-6
            theta_star_1 = theta_star_1_p * test(i, 3);
            theta_star_2 = theta_star_2_p * test(i, 2);
            break
        else
            if i == length(test)
                disp("Needs manual interference in delta_theta_star calc")
            end
        end

    end

    f_g_fdot_gdot = fg_out(R1, R2, delta_theta_star, p, mu);
    f = f_g_fdot_gdot(1);
    g = f_g_fdot_gdot(2);
    fdot = f_g_fdot_gdot(3);
    gdot = f_g_fdot_gdot(4);

    v_1_f = (R2 - f*R1)./g;
    v_2_i = fdot*R1 + gdot*v_1_f;

    delta_v_1 = v_1_f - V1;
    delta_v_2 = V2 - v_2_i;

    delta_v_tot = norm(abs(delta_v_1)) + norm(abs(delta_v_2));

    out = delta_v_tot;

end


function fx = lamberts_eqn_hyperbolic(a, s, c, TOF, mu, greater_than_180)
    alpha_h = 2 * asinh(s/(abs(a)));
    beta_h = 2 * asinh(sqrt((s-c)/(2*abs(a))));
    
    if greater_than_180
        fx = TOF - sqrt(a^3/mu) * (sinh(alpha_h) - alpha_h + (sinh(beta_h) - beta_h));
    else
        fx = TOF - sqrt(a^3/mu) * (sinh(alpha_h) - alpha_h - (sinh(beta_h) - beta_h));
    end

end


function fx = lamberts_eqn_elliptical(a, s, c, TOF, mu, TOF_min, greater_than_180)
    n = sqrt(mu/a^3);
    alpha_0 = 2 * asin(sqrt(s/(2*a)));
    beta_0 = 2 * asin(sqrt((s-c)/(2*a)));
    alpha_beta = lamberts_problem_alpha_beta_logic(TOF, TOF_min, greater_than_180, alpha_0, beta_0);
    alpha = alpha_beta(1);
    beta = alpha_beta(2);
    fx = TOF * n - (alpha - beta - (sin(alpha) - sin(beta)));
end
