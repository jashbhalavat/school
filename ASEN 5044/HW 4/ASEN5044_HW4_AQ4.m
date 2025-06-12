clear; clc; close all;

w = [0.1859, 0.0961, 0.1055, 0.2104, 0.0678, 0.1950, 0.1393];
a = [-4, -2, -0.75, 0, 3, 4, -0.9];
b = [-2, -1, 0, 1, 5.5, 6, 7];

x = -6:0.001:8;

for i = 1:length(x)
    px(i) = 0;
    for j = 1:length(w)
        px(i) = px(i) + w(j) * uni_dist_pdr(a(j), b(j), x(i));
    end
end

figure()
plot(x, px, 'LineWidth', 1.5)
xlabel("X values")
ylabel("Probability Density Funciton")
title("p(x)")

analytical_mean = 0;
analytical_var = 0;
for k = 1:length(w)
    analytical_mean = analytical_mean + w(k)*1/2*(a(k) + b(k));
    analytical_var = analytical_var + w(k) * (b(k) - a(k))^2/12;
end

numerical_mean = mean(px);
numerical_var = std(px)^2;

function uniform_pdf = uni_dist_pdr(a, b, x)
    if x >= a && x <= b
        uniform_pdf = 1/(b-a);
    else
        uniform_pdf = 0;
    end
end

%% Part c

e_x = 0;
var_x = 0;
for i = 1:length(px)-1
    e_x = x(i) * (x(i+1) - x(i)) * ((px(i+1) + px(i))/2) + e_x;
end

for i = 1:length(px)-1
    var_x = (x(i) - e_x)^2 * (x(i+1) - x(i)) * ((px(i+1) + px(i))/2) + var_x;
end

diff_ent = 0;
for i = 1:length(px)-1
    if px(i) == 0
        out = 0;
    else
        out = log(px(i));
    end
    diff_ent = -out * (x(i+1) - x(i)) * ((px(i+1) + px(i))/2) + diff_ent;
end


%%
% Part d
N1 = 100;
N2 = 1000;
N3 = 50000;
nbins = 1000;

var = 1:length(w);


C = randsample(var, 1, true, w);

N1_samples = part_d(N1, w, a, b, nbins);
N2_samples = part_d(N2, w, a, b, nbins);
N3_samples = part_d(N3, w, a, b, nbins);

function out = part_d(N, w, a, b, nbins)
    var = 1:length(w);
    for i = 1:N
        C = randsample(var, 1, true, w);
        samples(i) = a(C) + (b(C) - a(C))*rand;
    end

    figure()
    histogram(samples, nbins, 'Normalization','pdf')
    xlim([-6, 8])
    title("N = " + N)
    xlabel("X")
    ylabel("samples Count")

    out = samples;
end

%% Part e
e_x_n1 = 1/N1 * sum(N1_samples);
e_x_n2 = 1/N2 * sum(N2_samples);
e_x_n3 = 1/N3 * sum(N3_samples);

var_x_n1 = 1/N1 * sum((N1_samples - e_x_n1).^2);
var_x_n2 = 1/N2 * sum((N2_samples - e_x_n2).^2);
var_x_n3 = 1/N3 * sum((N3_samples - e_x_n3).^2);

h_x_n1 = 1/N1 * sum(-log(N1_samples));
h_x_n2 = 1/N2 * sum(-log(N2_samples));
h_x_n3 = 1/N3 * sum(-log(N3_samples));
