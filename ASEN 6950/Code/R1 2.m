function out = R1(ang)
    out = [1 0 0;
            0, cos(ang), sin(ang);
            0, -sin(ang), cos(ang)];
end