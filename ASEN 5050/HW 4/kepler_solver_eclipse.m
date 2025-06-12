function ecc_anomaly = kepler_solver_eclipse(time_past_periapsis, a, e, mu)
    % Numerically solve kepler's equation for eccentric anomaly using time
    % past periapsis, semi-major axis, eccentricity, system gravitational
    % parameter.
    % Following algorithm 2 in Vallado section 2.2.5

    tol = 1e-8;
    
    % Calculate period
    P = 2 * pi * sqrt(a^3/mu);
    
    % Mean anomaly - M = n(t-tp) = 2pi/p * (t-tp)
    M = (2*pi)/P * time_past_periapsis;
    
    % Using  series solution of E_N and truncating higher order
    % terms
    E_N = M + e*sin(M) + e^2/2 * sin(2*M);

    % E_N+1 can be obtained through the Newton-Raphson method
    E_N_plus_1 = E_N + (M - E_N + e*sin(E_N))/(1 - e*cos(E_N));

    % Keep applying the NR method until tolerance is met
    while abs(E_N_plus_1 - E_N) > tol
        E_N = E_N_plus_1;
        E_N_plus_1 = E_N + (M - E_N + e*sin(E_N))/(1 - e*cos(E_N));
    end

    % Output
    ecc_anomaly = E_N_plus_1;
end

