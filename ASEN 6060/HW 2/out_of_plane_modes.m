function out = out_of_plane_modes(mu, x_eq)
    % Calculate two out of plane modes for eq points
    lambda_pos = sqrt(u_zz(mu, x_eq));
    lambda_neg = -sqrt(u_zz(mu, x_eq));
    out = [lambda_pos, lambda_neg];
end