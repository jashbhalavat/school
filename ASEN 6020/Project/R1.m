function out = R1(theta)
    % DCM for a rotation around axis 1 by angle theta [rad]
    out = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
end