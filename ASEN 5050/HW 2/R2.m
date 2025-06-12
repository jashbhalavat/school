function out = R2(theta)
    % DCM for a rotation around axis 2 by angle theta [rad]
    out = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
end