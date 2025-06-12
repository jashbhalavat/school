function out = H(state0, statef, x_star_plus_delta)
    out = [statef(1) - state0(1);
            statef(2) - state0(2);
            statef(3) - state0(3);
            statef(4) - state0(4);
            statef(6) - state0(6);
            state0(2);
            state0(1) - x_star_plus_delta];
end