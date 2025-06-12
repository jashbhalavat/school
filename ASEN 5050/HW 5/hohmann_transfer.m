function out = hohmann_transfer(r_initial, r_final, mu)
    % This f'n outputs the total delta_v and time of flight for a hohmann
    % transfer
    % Algorithm 36 in Vallado

    a_trans = (r_initial + r_final)/2;

    v_initial = sqrt(mu/r_initial);
    v_final = sqrt(mu/r_final);

    v_trans_a = sqrt((2*mu)/(r_initial) - (mu)/(a_trans));
    v_trans_b = sqrt((2*mu)/(r_final) - (mu)/(a_trans));

    delta_v_a = v_trans_a - v_initial;
    delta_v_b = v_final - v_trans_b;
    
    delta_v_tot = abs(delta_v_a) + abs(delta_v_b);

    tof = pi * sqrt(a_trans^3 / mu);

    out = [delta_v_tot, tof];

end

