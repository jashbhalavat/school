function out = F(state0, statef)
    % Modified Constraint Vector
    out = statef(1:7) - state0(9:15); % x_1_f - x_2_0
end