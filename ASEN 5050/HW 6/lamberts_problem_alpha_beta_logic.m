function out = lamberts_problem_alpha_beta_logic(TOF, TOF_min, greater_than_180, alpha_0, beta_0)
    % Function that outputs appropriate alpha and beta depending on time of
    % flight (TOF) and choice of arc size (greater than 180 or no)
    % Inputs:
    %       TOF - Time of flight of transfer arc [s]
    %       TOF_min - Minimum energy arc TOF [s]
    %       greater_than_180 - bool that's true if chosen arc is >180 deg
    %       alpha_0, beta_0 - principal alpha, beta [rad]
    % Output:
    %       alpha, beta - actual alpha and beta [rad]

    % Logic from ASEN5050 Lambert's Problem alpha/beta table
    if TOF < TOF_min
        if greater_than_180
            % Purple arc
            alpha = alpha_0;
            beta = -beta_0;
        else
            % Blue arc
            alpha = alpha_0;
            beta = beta_0;
        end
    else
        if greater_than_180
            % Brown arc
            alpha = 2*pi - alpha_0;
            beta = -beta_0;
        else
            % Red arc
            alpha = 2*pi - alpha_0;
            beta = beta_0;
        end
    end

    out = [alpha, beta];
end