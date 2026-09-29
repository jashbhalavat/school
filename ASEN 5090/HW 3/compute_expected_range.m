function R1 = compute_expected_range(t_r, r_gps_tr, r_rx, c, ephem_data, prn, omega_E)
    % Computed expected range given Tr [WN, TOW], r_gps [ECEF, m], c [m/s]
    max_count = 5;

    % Initial measure of time taken for signal to get from satellite to
    % receiver
    dt = norm(r_gps_tr(1,:) - r_rx)/c;

    for j = 1:max_count
        for i = 1:length(t_r)
            % Step 2
            R(i) = norm(r_gps_tr(i,:) - r_rx);
    
            % Step 3
            t_t(i,1) = t_r(i,2) - R(i)/c;
        end
        t_t = [t_r(:,1), t_t];
    
        % Step 4
        [~, r_gps_tt] = eph2pvt2025(ephem_data, t_t, prn);
    
        % Step 5
        for i = 1:length(t_r)
            phi = omega_E * (t_r(i,2)-t_t(i,2));
            r_gps_tr(i,:) = ([cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1] * r_gps_tt(i,:)')';
    
            % Step 6
            R1(i,1) = norm(r_gps_tr(i,:) - r_rx);
        end

        % Calculate new signal travel time
        dt_new = R1(1)/c;

        % Compare signal travel times and break if within tolerance
        if abs(dt_new - dt) < 1e-12
            break
        end

        % If not within tolerance, update current dt
        dt = dt_new;
    end
end

