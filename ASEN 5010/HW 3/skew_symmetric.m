function out = skew_symmetric(vec)
    out = [0 -vec(3) vec(2);
            vec(3) 0 -vec(1);
            -vec(2) vec(1) 0];
end