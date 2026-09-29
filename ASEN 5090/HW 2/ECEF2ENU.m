function C_ECEF2ENU = ECEF2ENU(ref_lat_deg, ref_lon_deg)
    % Output transformation matrix that converts from ECEF to ENU.
    % Input angles are in deg

    % This matrix is from Slide 25 in Lecture 04
    C_ECEF2ENU = [-sind(ref_lon_deg), cosd(ref_lon_deg), 0;
                    -sind(ref_lat_deg)*cosd(ref_lon_deg), -sind(ref_lon_deg)*sind(ref_lat_deg), cosd(ref_lat_deg);
                    cosd(ref_lat_deg)*cosd(ref_lon_deg), cosd(ref_lat_deg)*sind(ref_lon_deg), sind(ref_lat_deg)];

end