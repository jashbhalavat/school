function out = sheppards_method(C)
    beta0_sq = 1/4 * (1 + trace(C));
    beta1_sq = 1/4 * (1 + 2*C(1,1) - trace(C));
    beta2_sq = 1/4 * (1 + 2*C(2,2) - trace(C));
    beta3_sq = 1/4 * (1 + 2*C(3,3) - trace(C));
    [val, max_idx] = max([beta0_sq, beta1_sq, beta2_sq, beta3_sq]);
    if max_idx == 1
        beta0 = sqrt(beta0_sq);
        beta1 = (C(2,3) - C(3,2))/(4*beta0);
        beta2 = (C(3,1) - C(1,3))/(4*beta0);
        beta3 = (C(1,2) - C(2,1))/(4*beta0);
    elseif max_idx == 2
        beta1 = sqrt(beta1_sq);
        beta0 = (C(2,3) - C(3,2))/(4*beta1);
        beta2 = (C(1,2) + C(2,1))/(4*beta1);
        beta3 = (C(3,1) + C(1,3))/(4*beta1);
    elseif max_idx == 3
        beta2 = sqrt(beta2_sq);
        beta0 = (C(3,1) - C(1,3))/(4*beta2);
        beta1 = (C(1,2) + C(2,1))/(4*beta2);
        beta3 = (C(2,3) + C(3,2))/(4*beta2);
    elseif max_idx == 4
        beta3 = sqrt(beta3_sq);
        beta0 = (C(1,2) - C(2,1))/(4*beta3);
        beta1 = (C(3,1) + C(1,3))/(4*beta3);
        beta2 = (C(2,3) + C(3,2))/(4*beta3);
    end
    out = [beta0, beta1, beta2, beta3];
    if beta0 < 0
        out = -out;
    end
end