% function out = F(state0, statef, x_1_des, x_4_des)
function out = F(state0, statef, x_2_des, init_mass)
    % Modified Constraint Vector
    F1 = statef(1:6) - state0(8:13); % x1f - x20
    F2 = [statef(8:13) - x_2_des; % x2f - x2des
            norm(state0(16:18))^2 - 1; % ||uhat2||^2 - 1
            state0(14) - init_mass]; % m_2_0 - init_mass
    % F2 = [statef(8:13) - state0(19:24); % x2f - x30
    %       statef(14) - state0(25); % m2f - m30
    %       norm(state0(16:18)) - 1]; % ||uhat2||^2 - 1
    % F3 = [statef(19:24) - state0(30:35); % x3f - x40
    %       norm(state0(27:29)) - 1]; % ||uhat3||^2 - 1
    % F4 = [statef(1:6) - x_1_des; % x10 - x1des
    %         statef(30:35) - x_4_des]; % x4f - x4des
    
    % out = [F1; F2; F3; F4];
    out = [F1; F2];
end