clear; clc; close all;

% ASEN 6020 - HW 1, Prob 4
% Jash Bhalavat
% Spring 2025


l = 10;
i_in = 90 * pi/180;

i1 = linspace(0, i_in, 1000);
i3 = linspace(0, i_in, 1000);

for j = 1:length(i1)
    for k = 1:length(i3)
        cost_dl(j,k) = dv_dl_norm(l, i1(j), i3(k), i_in);
    end
end
cost_be_range = dv_ib_norm(l, i_in);
% [min_cost_dl(p), index] = min(cost_dl(:));
% [min_i1_index(p), min_i3_index(p)] = ind2sub(size(cost_dl), index);

[min_cost_dl, index] = min(cost_dl(:));
[min_i1_index, min_i3_index] = ind2sub(size(cost_dl), index);

opt_i1 = i1(min_i1_index) * 180/pi;
opt_i3 = i1(min_i3_index) * 180/pi;

% avg_min_i1_index = mean(min_i1_index(2:end));
% avg_min_i3_index = mean(min_i3_index(2:end));
% 
% opt_i1 = i1(floor(avg_min_i1_index));
% opt_i3 = i3(floor(avg_min_i3_index));

% figure()
% plot(l_range, cost_be_range)
% hold on
% plot(l_range, min_cost_dl)
% hold off
% legend("BE", "DL")

figure()
contourf(i1, i3, cost_dl')
title("Contour plot of cost of transfer divided over 3 plane changes")
colorbar()
xlabel("\Delta i_1")
ylabel("\Delta i_3")

function out = dv_dl_norm(r, i1, i3, i_in)
    % di1 and di3 cost function
    dv1 = sqrt(1 + 2*r/(1+r) - 2*sqrt((2*r)/(1+r))*cos(i1));
    dv3 = sqrt(1 + 2*r/(1+r) - 2*sqrt((2*r)/(1+r))*cos(i3));
    dv2 = 2 * sqrt(2/(r*(1+r))) * sin((i_in - i1 - i3)/2);
    out = dv1 + dv2 + dv3;
end

function out = dv_ib_norm(r, i)
    out = 2 * (sqrt((2*r)/(1+r)) - 1 + sqrt((2)/(r*(1+r)))*sin(i/2));
end