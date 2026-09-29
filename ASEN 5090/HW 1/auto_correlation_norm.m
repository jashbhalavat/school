function out = auto_correlation_norm(n, x_k)
    % Normalized auto correlation function
    % n - shift
    % x_k - PRN associated with satellite k

    out = sum(x_k .* circshift(x_k, n))/1023;
    
end