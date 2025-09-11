function out = R3(theta)
    % DCM for a rotation around axis 3 by angle theta [rad]
    out = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
end