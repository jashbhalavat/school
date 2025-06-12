function out = euler_rotation_321(t1, t2, t3)
    % out = [cos(t2)*cos(t1), cos(t2)*sin(t1), -sin(t2);
    %         sin(t3)*sin(t2)*cos(t1) - cos(t3)*sin(t1), sin(t3)*sin(t2)*sin(t1) + cos(t3)*cos(t1), sin(t3)*cos(t2);
    %         cos(t3)*sin(t2)*cos(t1) + sin(t3)*sin(t1), cos(t3)*sin(t2)*sin(t1) - sin(t3)*cos(t1), cos(t3)*cos(t2)];
    out = R1(t3)*R2(t2)*R3(t1);
end