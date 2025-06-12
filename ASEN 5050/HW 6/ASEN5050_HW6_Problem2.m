clear; clc; close all;

% ASEN 5050 - HW 6
% Fall 2024, 10/31/2024
% Jash Bhalavat
% Problem 2

mu_sun = 1.32712428e11;
AU = 149597870.7;

% Earth to Mars
earth_ephem = readtable('HW6_Ephem_Earth.txt', 'ReadVariableNames', false, 'NumHeaderLines', 3);
for i = 1:height(earth_ephem)
    for j = 1:width(earth_ephem)
        earth_ephem_data(i, j) = str2double(earth_ephem{i, j});
    end
end

mars_ephem = readtable('HW6_Ephem_Mars.txt', 'ReadVariableNames', false, 'NumHeaderLines', 3);
for i = 1:height(mars_ephem)
    for j = 1:width(mars_ephem)
        mars_ephem_data(i, j) = str2double(mars_ephem{i, j});
    end
end

greater_than_180 = false;
delta_a = 0.1;

% Create for loops to go through every combination
for i = 1:length(earth_ephem_data)
    for j = 1:length(mars_ephem_data)
        earth_time = earth_ephem_data(i, 1);
        mars_time = mars_ephem_data(j, 1);

        % TOF in seconds
        TOF = (mars_time - earth_time) * 86400;

        % If TOF is positive, only then calculate the v_inf and store
        % departure and arrival times for plots later
        if TOF > 0
            R1 = [earth_ephem_data(i, 5), earth_ephem_data(i, 6), earth_ephem_data(i, 7)];
            R2 = [mars_ephem_data(j, 5), mars_ephem_data(j, 6), mars_ephem_data(j, 7)];
            V1 = [earth_ephem_data(i, 8), earth_ephem_data(i, 9), earth_ephem_data(i, 10)];
            V2 = [mars_ephem_data(j, 8), mars_ephem_data(j, 9), mars_ephem_data(j, 10)];
            v_inf = problem2_function(mu_sun, R1, V1, R2, V2, TOF, delta_a, greater_than_180);
            v_inf_earth(i, j) = v_inf(1);
            v_inf_mars(i, j) = v_inf(2);
            departure_time(i) = earth_time;
            arrival_time(j) = mars_time;
        else
            % If TOF is negative, disregard that combination by setting
            % v_inf and times to NaN
            v_inf_earth(i, j) = NaN;
            v_inf_mars(i, j) = NaN;
            departure_time(i) = NaN;
            arrival_time(j) = NaN;
        end
    end
end

%% Plot Porkchop Plot

% Set contour levels and plot
vec = 0:0.1:15;
figure()
contour(departure_time, arrival_time, v_inf_earth', vec)
colorbar
title("Porkchop Plot for v-infinity at Earth departure")
xlabel("Earth Departure Epoch [Julian Day TDB]")
ylabel("Mars Arrival Epoch [Julian Day TDB]")

figure()
contour(departure_time, arrival_time, v_inf_mars', vec)
colorbar
title("Porkchop Plot for v-infinity at Mars arrival")
xlabel("Earth Departure Epoch [Julian Day TDB]")
ylabel("Mars Arrival Epoch [Julian Day TDB]")

% Normalize time to start from first epoch to make it easier to read
% contour plots
departure_time_normalized = departure_time - departure_time(1);
arrival_time_normalized = arrival_time - arrival_time(1);
figure()
contour(departure_time_normalized, arrival_time_normalized, v_inf_earth', vec)
colorbar
title("Porkchop Plot for v-infinity at Earth departure")
xlabel("Earth Departure Epoch from 2005 June 20 at midnight [Julian Day TDB]")
ylabel("Mars Arrival Epoch from  2005 December 01 at midnight [Julian Day TDB]")

figure()
contour(departure_time_normalized, arrival_time_normalized, v_inf_mars', vec)
colorbar
title("Porkchop Plot for v-infinity at Mars arrival")
xlabel("Earth Departure Epoch from 2005 June 20 at midnight [Julian Day TDB]")
ylabel("Mars Arrival Epoch from  2005 December 01 at midnight [Julian Day TDB]")
