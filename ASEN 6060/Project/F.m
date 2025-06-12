function out = F(state0, statef, r_des_4)
    % Modified Constraint Vector
    out = [statef(1:3) - state0(8:10);
            statef(8:10) - state0(15:17);
            statef(15:17) - state0(22:24);
            statef(22:24) - r_des_4];
end