function out = compare_poincare_maps_beyond_moon(init, final, p2_pos_x)
    % Output two points on the poincare maps that are closest to each other
    % beyond the moon (towards the earth)
    % in 3 dimensions

    min = 1;
    for i = 1:length(init)
        for j = 1:length(final)
            if init(i,1) < p2_pos_x
                diff_norm = sqrt((init(i,1)-final(j,1))^2 + (init(i,2)-final(j,2))^2 + (init(i,3)-final(j,3))^2);
                if diff_norm < min
                    min = diff_norm;
                    i_min = i;
                    j_min = j;
                end
            end
        end
    end
    out(1) = min;
    out(2) = i_min;
    out(3) = j_min;
end