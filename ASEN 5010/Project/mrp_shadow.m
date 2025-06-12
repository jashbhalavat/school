function sigma_s = mrp_shadow(sigma)
    sigma_sq = dot(sigma, sigma);
    sigma_s = -sigma/sigma_sq;
end