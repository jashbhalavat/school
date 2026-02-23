function out = F(state0, statef, desired_stage_1, desired_stage_4)
    % Modified Constraint Vector
    out = [statef(1:6) - state0(8:13); % x_1_f - x_2_0
        statef(8:13) - state0(15:20); % x_2_f - x_3_0
        statef(15:20) - state0(22:27); % x_3_f - x_4_0
        statef(22:27) - desired_stage_4;
        state0(1:6) - desired_stage_1]; % x_4_f - x_f_des
end