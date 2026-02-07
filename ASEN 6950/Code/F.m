function out = F(state0, statef)
    % Modified Constraint Vector
    out = [statef(1) - state0(1); % x0
            statef(2) - state0(2); % y0
            statef(3) - state0(3); % z0
            statef(4) - state0(4); % x0_dot
            statef(6) - state0(6); % z0_dot
            state0(2)]; % y0
end