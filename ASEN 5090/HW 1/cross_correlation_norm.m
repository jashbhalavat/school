function out = cross_correlation_norm(n, x_k, x_l)
    % Normalized cross correlation function
    % n - shift
    % x_k - PRN associated with satellite k
    % x_l - PRN associated with satellite l

    out = sum(x_k .* circshift(x_l, [0, n]))/1023;
    
end