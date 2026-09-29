function [AZ, EL, RANGE] = compute_azelrange(userECEF, satECEF)
    % Computer azimuth [deg], elevation [deg], range [m] of GPS satellite
    % Inputs:
    %   userECEF - user ECEF position [1x3], [meters]
    %   satECEF - satellite ECEF position [nx3], [meters]
    % Outputs:
    %   AZ - Azimuth [nx1], [deg]
    %   EL - Elevation [nx1], [deg]
    %   RANGE - Range [nx1], [meters]

    % Convert user ecef to lla and use that to calculate C_ECEF2ENU
    % transformation matrix
    user_lla = ecef2lla(userECEF);
    C_ECEF2ENU_user = ECEF2ENU(user_lla(1), user_lla(2));

    % number of position vectors
    n = size(satECEF,1);

    % Create output vectors
    AZ = NaN(n,1);
    EL = NaN(n,1);
    RANGE = NaN(n,1);

    for i = 1:n
        % calculate user to satellite vector in ECEF and convert to ENU
        r_user2sat_ecef = satECEF(i,:) - userECEF;
        r_user2sat_enu = (C_ECEF2ENU_user * r_user2sat_ecef')';

        % Range is just magnitude of the vector
        RANGE(i) = norm(r_user2sat_enu);

        % Azimuth and elevation equations from slide 34, lecture 05
        AZ(i) = mod(atan2d(r_user2sat_enu(1), r_user2sat_enu(2)), 360);
        EL(i) = asind(r_user2sat_enu(3)/RANGE(i));
    end
end