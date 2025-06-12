function out = bi_elliptic_transfer(r_initial, r_final, r_b, mu)
    % This f'n outputs the total delta_v and time of flight for a
    % bi-elliptic transfer
    % Algorithm 37 in Vallado

    a_trans_1 = (r_initial + r_b)/2;
    a_trans_2 = (r_b + r_final)/2;

    v_initial = sqrt(mu/r_initial);
    v_trans_1_a = sqrt((2*mu)/(r_initial) - (mu)/(a_trans_1));
    v_trans_1_b = sqrt((2*mu)/(r_b) - (mu)/(a_trans_1));
    v_trans_2_b = sqrt((2*mu)/(r_b) - (mu)/(a_trans_2));
    v_trans_2_c = sqrt((2*mu)/(r_final) - (mu)/(a_trans_2));
    v_final = sqrt(mu/r_final);

    delta_v_a = v_trans_1_a - v_initial;
    delta_v_b = v_trans_2_b - v_trans_1_b;
    delta_v_c = v_final - v_trans_2_c;
    
    delta_v_tot = abs(delta_v_a) + abs(delta_v_b) + abs(delta_v_c);

    tof = pi * sqrt(a_trans_1^3 / mu) + pi * sqrt(a_trans_2^3 / mu);

    out = [delta_v_tot, tof];

end

