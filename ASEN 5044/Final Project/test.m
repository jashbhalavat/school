clear; clc; close all;

% ASEN 5044
% Final Project

% Given constants
rE = 6378; % km
omegaE = 2*pi/86400; % rad/s
alt = 300; % km
r0 = rE + alt;
mu = 398600; % km^3/s^2
stations = 12;

% Initial conditions
perturb_x0 = [0;0.075;0;-0.021];
x0 = [r0; 0; 0; r0*sqrt(mu/r0^3)] + perturb_x0;

dt = 10; % s
t0 = 0;
tp = 2*pi*sqrt(r0^3/mu);
tf = 1400*dt;
tspan = t0:dt:tf;

%% Linearized Dynamics

tspan = t0:dt:tf;

dx = NaN(4,length(tspan));
xnom = NaN(4,length(tspan));
xtot = NaN(4,length(tspan));
A = NaN(4,4,length(tspan));
F = NaN(4,4,length(tspan));

dx(:,1) = perturb_x0;
xnom(:,1) = x0;
xtot(:,1) = xnom(:,1) + dx(:,1);

A(:,:,1) = [0 1 0 0; -mu*(xnom(1,1)^2+xnom(3,1)^2)^(-3/2) + 3*mu*xnom(1,1)^2*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0 3*mu*xnom(1,1)*xnom(3,1)*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0; ...
        0 0 0 1; 3*mu*xnom(1,1)*xnom(3,1)*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0 -mu*(xnom(1,1)^2+xnom(3,1)^2)^(-3/2) + 3*mu*xnom(3,1)^2*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0];
F(:,:,1) = eye(4) + A(:,:,1)*dt;

for i=2:length(tspan)
    % Nominal state
    x = r0*cos(2*pi/tp*tspan(i));
    y = r0*sin(2*pi/tp*tspan(i));
    vx = -2*pi/tp*y;
    vy = 2*pi/tp*x;
    xnom(:,i) = [x; vx; y; vy];

    % perturbations linearized about nominal point
    A(:,:,i) = [0 1 0 0; -mu*(xnom(1,i)^2+xnom(3,i)^2)^(-3/2) + 3*mu*xnom(1,i)^2*(xnom(1,i)^2+xnom(3,i)^2)^(-5/2) 0 3*mu*xnom(1,i)*xnom(3,i)*(xnom(1,i)^2+xnom(3,i)^2)^(-5/2) 0; ...
        0 0 0 1; 3*mu*xnom(1,i)*xnom(3,i)*(xnom(1,i)^2+xnom(3,i)^2)^(-5/2) 0 -mu*(xnom(1,i)^2+xnom(3,i)^2)^(-3/2) + 3*mu*xnom(3,i)^2*(xnom(1,i)^2+xnom(3,i)^2)^(-5/2) 0];
    F(:,:,i) = eye(4) + dt * A(:,:,i);
    dx(:,i) = F(:,:,i-1)*dx(:,i-1);
    xtot(:,i) = xnom(:,i) + dx(:,i);
end

%% NL ODE45

% ode45 solution
tol = 1e-12;
options = odeset('Stats','off','RelTol',tol,'AbsTol',tol);
[tout, xout] = ode45(@(tout,xout) eom(xout,mu),tspan,x0,options);

function xdot = eom(x,mu)
    r = sqrt(x(1)^2+x(3)^2);
    xdot = [x(2); -mu*x(1)/r^3; x(4); -mu*x(3)/r^3];
end

% tracking station positions using ode45
Xi = NaN(length(tspan),stations);
Yi = NaN(length(tspan),stations);
Xdoti = NaN(length(tspan),stations);
Ydoti = NaN(length(tspan),stations);
thetai = NaN(length(tspan),stations);
phii_nom = NaN(length(tspan),stations);
rhoi_nom = NaN(length(tspan),stations);
rhoidot_nom = NaN(length(tspan),stations);
satID_nom = NaN(length(tspan),stations);

for i=1:stations
   thetai(1,i) = (i-1)*pi/6;
   Xi(1,i) = rE*cos(thetai(1,i));
   Yi(1,i) = rE*sin(thetai(1,i));
   Xdoti(1,i) = -omegaE*Yi(1,i);
   Ydoti(1,i) = omegaE*Xi(1,i);
   phii_nom(1,i) = atan2(xout(1,3)-Yi(1,i), xout(1,1)-Xi(1,i));
   rhoi_nom(1,i) = sqrt((xout(1,1)-Xi(1,i))^2 + (xout(1,3)-Yi(1,i))^2);
   rhoidot_nom(1,i) = ((xout(1,1)-Xi(1,i))*(xout(1,2)-Xdoti(1,i)) + ((xout(1,3)-Yi(1,i))*(xout(1,4)-Ydoti(1,i))))/rhoi_nom(1,i);
   
   sc_from_gs = [xnom(1,1) - Xi(1,i), xnom(3,1) - Yi(1,i)];
   angle = ang_diff(sc_from_gs, [Xi(1,i), Yi(1,i)]);
   
   philower = wrapToPi(-pi/2+thetai(1,i));
   phiupper = wrapToPi(pi/2+thetai(1,i));
   if philower < phiupper
       if ~(phii_nom(1,i) >= philower && phii_nom(1,i) <= phiupper)
           satID_nom(1,i) = NaN;
           rhoi_nom(1,i) = NaN;
           rhoidot_nom(1,i) = NaN;
           phii_nom(1,i) = NaN;
       end
   else
       if ~(phii_nom(1,i) >= philower || phii_nom(1,i) <= phiupper)
           satID_nom(1,i) = NaN;
           rhoi_nom(1,i) = NaN;
           rhoidot_nom(1,i) = NaN;
           phii_nom(1,i) = NaN;
       end
   end
end

for k=2:length(tspan)
    for i=1:stations
        satID_nom(k,i) = i;
        Xi(k,i) = rE*cos(omegaE*tspan(k)+thetai(1,i));
        Yi(k,i) = rE*sin(omegaE*tspan(k)+thetai(1,i));
        Xdoti(k,i) = -omegaE*Yi(k,i);
        Ydoti(k,i) = omegaE*Xi(k,i);
        thetai(k,i) = atan2(Yi(k,i),Xi(k,i));
        phii_nom(k,i) = atan2(xout(k,3)-Yi(k,i), xout(k,1)-Xi(k,i));
        rhoi_nom(k,i) = sqrt((xout(k,1)-Xi(k,i))^2 + (xout(k,3)-Yi(k,i))^2);
        rhoidot_nom(k,i) = ((xout(k,1)-Xi(k,i))*(xout(k,2)-Xdoti(k,i)) + ((xout(k,3)-Yi(k,i))*(xout(k,4)-Ydoti(k,i))))/rhoi_nom(k,i);
        
        philower = wrapToPi(-pi/2+thetai(k,i));
        phiupper = wrapToPi(pi/2+thetai(k,i));
        if philower < phiupper
            if ~(phii_nom(k,i) >= philower && phii_nom(k,i) <= phiupper)
                satID_nom(k,i) = NaN;
                rhoi_nom(k,i) = NaN;
                rhoidot_nom(k,i) = NaN;
                phii_nom(k,i) = NaN;
            end
        else
            if ~(phii_nom(k,i) >= philower || phii_nom(k,i) <= phiupper)
                satID_nom(k,i) = NaN;
                rhoi_nom(k,i) = NaN;
                rhoidot_nom(k,i) = NaN;
                phii_nom(k,i) = NaN;
            end
        end
    end
end



%% linearized solution

thetai = NaN(length(tspan),stations);
phii = NaN(length(tspan),stations);
rhoi = NaN(length(tspan),stations);
rhoidot = NaN(length(tspan),stations);
satID = NaN(length(tspan),stations);
ynom = NaN(3,stations, length(tspan));
ytot = NaN(3,stations, length(tspan));
dy = NaN(3,stations, length(tspan));
xi = NaN(4, stations, length(tspan));

for k = 1:length(tspan)
    for i=1:stations
       satID(k,i) = i;
       thetai(1,i) = (i-1)*pi/6;
       Xi = rE*cos(omegaE*tspan(k)+thetai(1,i));
       Yi = rE*sin(omegaE*tspan(k)+thetai(1,i));
       Xdoti = -omegaE*Yi;
       Ydoti = omegaE*Xi;
       xi(:,i,k) = [Xi; Xdoti; Yi; Ydoti];
    
       ynom(:,i,k) = h(xnom(:,k), xi(:,i,k));
       dy(:,i,k) = y_perturb(ynom(:,i,k), xi(:,i,k), xnom(:,k), dx(:,k));
       ytot(:,i,k) = ynom(:,i,k) + dy(:,i,k);
       rhoi(k,i) = ytot(1,i,k);
       rhoidot(k,i) = ytot(2,i,k);
       phii(k,i) = ytot(3,i,k);
       
       if k > 1
            thetai(k,i) = atan2(Yi,Xi);
       end
        
       if isnan(satID_nom(k, i))
            satID(k,i) = NaN;
            rhoi(k,i) = NaN;
            rhoidot(k,i) = NaN;
            phii(k,i) = NaN;
       end
   end
end


function out = h(x_star, x_i)
    out(1,1) = sqrt((x_star(1)-x_i(1))^2 + (x_star(3)-x_i(3))^2);
    out(2,1) = ((x_star(1)-x_i(1))*(x_star(2)-x_i(2)) + (x_star(3)-x_i(3))*(x_star(4)-x_i(4)))/out(1);
    out(3,1) = atan2((x_star(3)-x_i(3)),(x_star(1)-x_i(1)));
end

function dy = y_perturb(y_nom, X_gs, xnom, dx)
    xd = xnom - X_gs; % X_diff
    rho_i = y_nom(1);
    rho_dot_i = y_nom(2);
    rho_dot_i_num = xd(1)*xd(2) + xd(3)*xd(4);
    dh1_dx1 = xd(1)/rho_i;
    dh1_dx3 = xd(3)/rho_i;
    dh2_dx1 = xd(2)/rho_i - rho_dot_i*dh1_dx1/rho_i;
    dh2_dx3 = xd(4)/rho_i - rho_dot_i*dh1_dx3/rho_i;
    dh2_dx2 = xd(1)/rho_i;
    dh2_dx4 = xd(3)/rho_i;
    dh3_dx1 = 1/(1 + (xd(3)/xd(1))^2) * -(xd(3)/xd(1)^2);
    dh3_dx3 = 1/(1 + (xd(3)/xd(1))^2) * (1/xd(1));
    dh_dx = [dh1_dx1, 0, dh1_dx3, 0;
        dh2_dx1, dh2_dx2, dh2_dx3, dh2_dx4;
        dh3_dx1, 0, dh3_dx3, 0];
    dy = dh_dx * dx;
end

function angle = ang_diff(V1, V2)
    angle = acos(dot(V1, V2) / (norm(V1) * norm(V2)));
end


% Plot perturbations
figure()
subplot(4,1,1)
plot(tspan,dx(1,:))
title('Linearized perturbations')
xlabel('Time (s)')
ylabel('\deltax_1 (m)')
subplot(4,1,2)
plot(tspan,dx(2,:))
xlabel('Time (s)')
ylabel('\deltax_2 (m/s)')
subplot(4,1,3)
plot(tspan,dx(3,:))
xlabel('Time (s)')
ylabel('\deltax_3 (m)')
subplot(4,1,4)
plot(tspan,dx(4,:))
xlabel('Time (s)')
ylabel('\deltax_4 (m/s)')

% Plot full state
figure()
subplot(4,1,1)
plot(tspan,xtot(1,:))
title('Linearized solution')
xlabel('Time (s)')
ylabel('x_1 (m)')
subplot(4,1,2)
plot(tspan,xtot(2,:))
xlabel('Time (s)')
ylabel('x_2 (m/s)')
subplot(4,1,3)
plot(tspan,xtot(3,:))
xlabel('Time (s)')
ylabel('x_3 (m)')
subplot(4,1,4)
plot(tspan,xtot(4,:))
xlabel('Time (s)')
ylabel('x_4 (m/s)')

% Plot ground stations
figure()
subplot(4,1,1)
plot(tspan,rhoi, 'Marker','x')
title('Linearized ground stations')
xlabel('Time (s)')
ylabel('\rho^i (m)')
subplot(4,1,2)
plot(tspan,rhoidot, 'Marker','o')
xlabel('Time (s)')
ylabel('\rho^idot (m/s)')
subplot(4,1,3)
plot(tspan,phii,'+')
xlabel('Time (s)')
ylabel('\phi^i (radians)')
subplot(4,1,4)
plot(tspan,satID,'+')
xlabel('Time (s)')
ylabel('Station ID')