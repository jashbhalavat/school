function out = F(state0, statef, x_1_des, x_4_des)
    % Modified Constraint Vector
    % F1 = [state0(1:6) - x_1_des; % x10 - x1des 
    %         statef(1:6) - state0(9:14)]; % x1f - x20
    F1 = [state0(1:6) - x_1_des; % x_1_0 - x_1_des
            state0(7) - 1; % m_1_0 - m_1_des
            statef(1:6) - state0(9:14); % x1f - x20
            statef(7) - state0(15)]; % m_1_f - m_2_0
    F2 = [statef(9:14) - state0(20:25); % x_2_f - x_3_0
            statef(15) - state0(26); % m_2_f - m_3_0
            norm(state0(16:18))^2 - 1]; % ||uhat2||^2 - 1
    F3 = [statef(20:25) - state0(31:36); % x_3_f - x_4_0
            statef(26) - state0(37); % m_3_f - m_4_0
            norm(state0(27:29))^2 - 1]; % ||uhat3||^2 - 1
    F4 = statef(31:36) - x_4_des; % x_4_f - x_4_des

    out = [F1; F2; F3; F4];
end