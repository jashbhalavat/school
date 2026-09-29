function XG = ca_code(final_epoch, S1, S2)
    % Generate C/A code for GPS satellites.
    % final_epoch - Code will be generated for 1:final_epoch
    % S1, S2 - Phase selections. Can be found - https://archive.gps.gov/technical/icwg/IS-GPS-200D.pdf

    G1 = ones(1,10);
    G2 = ones(1,10);

    for i = 1:final_epoch
        G1_newbit = xor(G1(3), G1(10));
        G2_newbit = xor(G2(2), xor(G2(3), xor(G2(6), xor(G2(8), xor(G2(9), G2(10))))));
        G2i = xor(G2(S1), G2(S2));
        XG(i) = xor(G2i, G1(10));
    
        G1 = [G1_newbit G1(1:9)];
        G2 = [G2_newbit G2(1:9)];
    end

end