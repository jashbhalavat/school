function vec = cell2vec(input_cell)
    % Convert cell to a 1d column vector
    vec = [];

    for i = 1:length(input_cell)
        vec = [vec; input_cell{i}];
    end
end