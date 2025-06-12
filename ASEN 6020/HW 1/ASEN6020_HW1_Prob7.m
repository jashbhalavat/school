close all; clear; clc;

% ASEN 6020, HW 1, Problem 7
% Jash Bhalavat
% Spring 2025

v_inf_range = linspace(0, 2, 100);
R_range = linspace(0, 10, 1000);

for i = 1:length(v_inf_range)
    for j = 1:length(R_range)
        % dva = v_inf_range(i)-1;
        dva = (sqrt(2)-1);
        R = R_range(j);
        if R < 1
            dvb = 1 - sqrt((2*R)/(1+R)) + v_inf_range(i) - sqrt((2)/(R*(1+R)));
            % dvb = 1 - sqrt((2*R)/(1+R)) + sqrt(2/R) - sqrt((2)/(R*(1+R)));
            dvc = 10000;
        elseif R > 1
            dvc = - 1 + sqrt((2*R)/(1+R)) + v_inf_range(i) - sqrt((2)/(R*(1+R)));
            % dvc = - 1 + sqrt((2*R)/(1+R)) + sqrt(2/R) - sqrt((2)/(R*(1+R)));
            dvb = 10000;
        end
        costs = [dva, dvb, dvc];
        [m, ind] = min(costs);
        % if ind == 3
        %     disp("test")
        % end
        cost(i,j) = ind;
    end
end

figure()
contourf(v_inf_range, R_range, cost')



