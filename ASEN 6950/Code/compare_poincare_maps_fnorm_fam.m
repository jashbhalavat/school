function out = compare_poincare_maps_fnorm_fam(init, final)
    % Output the 25 pairs of points on the poincare maps that are closest 
    % to each other using constraint vector norm.
    %
    % Output format (matrix of size 25 x 3):
    % Column 1: Difference Norm Value
    % Column 2: 'init' row index (i)
    % Column 3: 'final' row index (j)

    len_init = size(init, 1);
    len_final = size(final, 1);
    total_combinations = len_init * len_final;

    % Preallocate arrays for speed
    all_norms = zeros(total_combinations, 1);
    all_i = zeros(total_combinations, 1);
    all_j = zeros(total_combinations, 1);

    idx = 1;
    for i = 1:len_init
        for j = 1:len_final
            % Compute the difference vector
            diff_vec = [init(i,1)-final(j,1), ...
                        init(i,2)-final(j,2), ...
                        init(i,3)-final(j,3), ...
                        init(i,end)-final(j,end-1)];
                    
            all_norms(idx) = norm(diff_vec);
            all_i(idx) = i;
            all_j(idx) = j;
            idx = idx + 1;
        end
    end

    % Sort all results by the norm value in ascending order
    [sorted_norms, sort_idx] = sort(all_norms, 'ascend');

    % Determine how many items to return (exactly 25, or fewer if matrix is small)
    num_outputs = min(25, total_combinations);

    % Extract the top elements using the sorted index tracker
    out = zeros(num_outputs, 3);
    out(:, 1) = sorted_norms(1:num_outputs);
    out(:, 2) = all_i(sort_idx(1:num_outputs));
    out(:, 3) = all_j(sort_idx(1:num_outputs));
end
