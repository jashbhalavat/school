function print_array(array)
    fprintf("%s - [", inputname(1))
    for i = 1:length(array)
        if i == length(array)
            if abs(array(i)) < 1e-12
                fprintf("0.0]^T\n")
            else
                fprintf("%0.14f]^T\n", array(i));
            end
        else
            if abs(array(i)) < 1e-12
                fprintf("0.0, ")
            else
                fprintf("%0.14f, ", array(i));
            end
        end
    end
end