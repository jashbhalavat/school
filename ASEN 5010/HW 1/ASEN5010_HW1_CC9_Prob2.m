clear; clc;

D2R = pi/180;

t1 = 40*D2R;
t2 = 30*D2R;
t3 = 80*D2R;

y0 = [t1; t2; t3];

tspan = [0:60]; % sec

[t, y] = ode45(@diff_eq, tspan, y0);

norm_at_42 = sqrt(y(43,1)^2 + y(43,2)^2 + y(43,3)^2);

function out = omega_B(t)
    out = 20 * [sin(0.1*t); 0.01; cos(0.1*t)] * pi/180;
end

function out = mult_mat(y)
    out = 1/cos(y(2)) * [0, sin(y(3)), cos(y(3)); 
        0, cos(y(3))*cos(y(2)), -sin(y(3))*cos(y(2)); 
        cos(y(2)), sin(y(3))*sin(y(2)), cos(y(3))*sin(y(2))];
end

function ydot = diff_eq(t, y)
    ydot = mult_mat(y) * omega_B(t);
end


