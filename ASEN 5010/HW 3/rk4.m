function yout = rk4(ode_fun, tspan, y0)

    % a = [0; 1/2; 1/2; 1];
    % b = [0 0 0;
    %     1/2 0 0;
    %     0 1/2 0;
    %     0 0 1];
    % c = [1/6; 1/3; 1/3; 1/6];

    yout = y0;

    h = tspan(2) - tspan(1);

    for i = 2:length(tspan)
        t = tspan(i);
        y = yout(i-1,:);
        
        f1 = ode_fun(t, y);
        f2 = ode_fun(t + 1/2*h, y + 1/2*h*f1);
        f3 = ode_fun(t + 1/2*h, y + 1/2*h*f2);
        f4 = ode_fun(t + h, y + h*f3);
        yout(i,:) = y + h * (1/6*f1 + 1/3*f2 + 1/3*f3 + 1/6*f4)';

        mrp_norm = norm(yout(i,:));
        if mrp_norm > 1
            yout(i,:) = -yout(i,:)./mrp_norm^2;
        end
end