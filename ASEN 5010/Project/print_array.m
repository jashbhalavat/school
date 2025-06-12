function print_array(array)
    if size(array, 1) == 1 || size(array, 2) == 1
        fprintf("%s - ", inputname(1))
        for i = 1:length(array)
            fprintf("%0.8f ", array(i))
        end
        fprintf("\n");
    else
        fprintf("%s - ", inputname(1))
        for i = 1:size(array, 1)
            for j = 1:size(array, 2)
                fprintf("%0.8f ", array(i, j))
        end
        fprintf("\n");
    end
end
