clear; clc; close all;

% Slick's offer - z1 (thousands of $)
% Prad Bitt's offer - z2 (thousands of $)

z0 = [100; 20];
z1 = [43.6658; 39.2815];
z2 = [40.5785; 40.3382];
z3 = [40.4093; 40.3961];
z4 = [40.4; 40.3993];
z5 = [40.3995; 40.3995];

z = [z0, z1, z2, z3, z4, z5];

for i = 2:6
    lambda(i) = (z(1,i) - z(1,i-1))/(z(1,i-1) - z(2,i-1));
    mu(i) = (z(2,i) - z(2,i-1))/(z(1,i-1) - z(2,i-1));
end

for i = 1:length(z)-1
    y(:,i) = [z(1,i+1) - z(1,i); z(2,i+1) - z(2,i)];
end
y = reshape(y, [10, 1]);

H = [];
for i = 1:length(z)-1
    H_temp = [z(1,i) - z(2,i), 0; 0, z(1,i)-z(2,i)];
    H = [H; H_temp];
end

gram = H' * H;
x = inv(gram) * H' * y;
