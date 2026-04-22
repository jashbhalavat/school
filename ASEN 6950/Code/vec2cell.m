function cell_out = vec2cell(V_vec, V_config)
    % Convert V in 1d col vector to cell form

    cell_out = cell(length(V_config),1);
    
    start_index = 1;
    for i = 1:length(V_config)
        if V_config(i) == 0
            cell_out{i} = V_vec(start_index:start_index+7);
            start_index = start_index + 8;
        elseif V_config(i) == 1
            cell_out{i} = V_vec(start_index:start_index+10);
            start_index = start_index + 11;
        end
    end
end