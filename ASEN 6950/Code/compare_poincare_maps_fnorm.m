function out = compare_poincare_maps_fnorm(init, final)
    % Output two points on the poincare maps that are closest to each other
    % using constraint vector norm

    min = 1;
    for i = 1:length(init)
        for j = 1:length(final)
            diff_vec = [init(i,1)-final(j,1), init(i,2)-final(j,2), init(i,3)-final(j,3), init(i,end)-final(j,end)];
            diff_norm = norm(diff_vec);
            if diff_norm < min
                min = diff_norm;
                i_min = i;
                j_min = j;
            end
        end
    end
    out(1) = min;
    out(2) = i_min;
    out(3) = j_min;
end