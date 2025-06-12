clear
clc

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
x0 = [r0; 0; 0; r0*sqrt(mu/r0^3)];

dt = 10; % s
t0 = 0;
tp = 2*pi*sqrt(r0^3/mu);
tf = 1400*dt;

tspan = t0:dt:tf;

% ode45 solution
tol = 1e-8;
options = odeset('Stats','off','RelTol',tol,'AbsTol',tol);
[tout, xout] = ode45(@(tout,xout) eom(xout,mu),tspan,x0,options);

function xdot = eom(x,mu)
    r = sqrt(x(1)^2+x(3)^2);
    xdot = [x(2); -mu*x(1)/r^3; x(4); -mu*x(3)/r^3];
end

% ode45 plots
figure()
subplot(4,1,1)
plot(tout,xout(:,1))
title('ode45 Simulation')
xlabel('Time (s)')
ylabel('x_1 (m)')
subplot(4,1,2)
plot(tout,xout(:,2))
xlabel('Time (s)')
ylabel('x_2 (m/s)')
subplot(4,1,3)
plot(tout,xout(:,3))
xlabel('Time (s)')
ylabel('x_3 (m)')
subplot(4,1,4)
plot(tout,xout(:,4))
xlabel('Time (s)')
ylabel('x_4 (m/s)')


% tracking station positions using ode45
Xi = NaN(length(tout),stations);
Yi = NaN(length(tout),stations);
Xdoti = NaN(length(tout),stations);
Ydoti = NaN(length(tout),stations);
thetai = NaN(length(tout),stations);
phii = NaN(length(tout),stations);
rhoi = NaN(length(tout),stations);
rhoidot = NaN(length(tout),stations);
satID = NaN(length(tout),stations);

for i=1:stations
   thetai(1,i) = (i-1)*pi/6;
   Xi(1,i) = rE*cos(thetai(1,i));
   Yi(1,i) = rE*sin(thetai(1,i));
   Xdoti(1,i) = -omegaE*Yi(1,i);
   Ydoti(1,i) = omegaE*Xi(1,i);
   phii(1,i) = atan2(xout(1,3)-Yi(1,i), xout(1,1)-Xi(1,i));
   rhoi(1,i) = sqrt((xout(1,1)-Xi(1,i))^2 + (xout(1,3)-Yi(1,i))^2);
   rhoidot(1,i) = ((xout(1,1)-Xi(1,i))*(xout(1,2)-Xdoti(1,i)) + ((xout(1,3)-Yi(1,i))*(xout(1,4)-Ydoti(1,i))))/rhoi(1,i);
   philower = wrapToPi(-pi/2+thetai(1,i));
   phiupper = wrapToPi(pi/2+thetai(1,i));
   if philower < phiupper
       if ~(phii(1,i) >= philower && phii(1,i) <= phiupper)
           satID(1,i) = NaN;
           rhoi(1,i) = NaN;
           rhoidot(1,i) = NaN;
           phii(1,i) = NaN;
       end
   else
       if ~(phii(1,i) >= philower || phii(1,i) <= phiupper)
           satID(1,i) = NaN;
           rhoi(1,i) = NaN;
           rhoidot(1,i) = NaN;
           phii(1,i) = NaN;
       end
   end
end

for k=2:length(tout)
    for i=1:stations
        satID(k,i) = i;
        Xi(k,i) = rE*cos(omegaE*tout(k)+thetai(1,i));
        Yi(k,i) = rE*sin(omegaE*tout(k)+thetai(1,i));
        Xdoti(k,i) = -omegaE*Yi(k,i);
        Ydoti(k,i) = omegaE*Xi(k,i);
        thetai(k,i) = atan2(Yi(k,i),Xi(k,i));
        phii(k,i) = atan2(xout(k,3)-Yi(k,i), xout(k,1)-Xi(k,i));
        rhoi(k,i) = sqrt((xout(k,1)-Xi(k,i))^2 + (xout(k,3)-Yi(k,i))^2);
        rhoidot(k,i) = ((xout(k,1)-Xi(k,i))*(xout(k,2)-Xdoti(k,i)) + ((xout(k,3)-Yi(k,i))*(xout(k,4)-Ydoti(k,i))))/rhoi(k,i);
        philower = wrapToPi(-pi/2+thetai(k,i));
        phiupper = wrapToPi(pi/2+thetai(k,i));
        if philower < phiupper
            if ~(phii(k,i) >= philower && phii(k,i) <= phiupper)
                satID(k,i) = NaN;
                rhoi(k,i) = NaN;
                rhoidot(k,i) = NaN;
                phii(k,i) = NaN;
            end
        else
            if ~(phii(k,i) >= philower || phii(k,i) <= phiupper)
                satID(k,i) = NaN;
                rhoi(k,i) = NaN;
                rhoidot(k,i) = NaN;
                phii(k,i) = NaN;
            end
        end
    end
end

figure()
subplot(4,1,1)
plot(tout,rhoi)
title('Tracking stations with ode45')
xlabel('Time (s)')
ylabel('\rho^i (m)')
subplot(4,1,2)
plot(tout,rhoidot)
xlabel('Time (s)')
ylabel('\rho^idot (m/s)')
subplot(4,1,3)
plot(tout,phii,'+')
xlabel('Time (s)')
ylabel('\phi^i (radians)')
subplot(4,1,4)
plot(tout,satID,'+')
xlabel('Time (s)')
ylabel('Station ID')

% linearized solution

tspan = t0:dt:tf;

dx = NaN(4,length(tspan));
xnom = NaN(4,length(tspan));
xtot = NaN(4,length(tspan));
A = NaN(4,4,length(tspan));
F = NaN(4,4,length(tspan));

dx(:,1) = perturb_x0;
% dx(:,1) = [0;0;0;0];
xnom(:,1) = x0;
xtot(:,1) = xnom(:,1) + dx(:,1);

A(:,:,1) = [0 1 0 0; -mu*(xnom(1,1)^2+xnom(3,1)^2)^(-3/2) + 3*mu*xnom(1,1)^2*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0 3*mu*xnom(1,1)*xnom(3,1)*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0; ...
        0 0 0 1; 3*mu*xnom(1,1)*xnom(3,1)*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0 -mu*(xnom(1,1)^2+xnom(3,1)^2)^(-3/2) + 3*mu*xnom(3,1)^2*(xnom(1,1)^2+xnom(3,1)^2)^(-5/2) 0];
F(:,:,1) = expm(A(:,:,1)*dt);

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
    F(:,:,i) = expm(A(:,:,i)*dt);
    dx(:,i) = F(:,:,i)*dx(:,i-1);
    xtot(:,i) = xnom(:,i) + dx(:,i);
end

% ground station positions using linearized model
for k=2:length(tspan)
    for i=1:stations
        satID(k,i) = i;
        Xi(k,i) = rE*cos(omegaE*tspan(k)+thetai(1,i));
        Yi(k,i) = rE*sin(omegaE*tspan(k)+thetai(1,i));
        Xdoti(k,i) = -omegaE*Yi(k,i);
        Ydoti(k,i) = omegaE*Xi(k,i);
        thetai(k,i) = atan2(Yi(k,i),Xi(k,i));
        phii(k,i) = atan2(xtot(3,k)-Yi(k,i), xtot(1,k)-Xi(k,i));
        rhoi(k,i) = sqrt((xtot(1,k)-Xi(k,i))^2 + (xtot(3,k)-Yi(k,i))^2);
        rhoidot(k,i) = ((xtot(1,k)-Xi(k,i))*(xtot(2,k)-Xdoti(k,i)) + ((xtot(3,k)-Yi(k,i))*(xtot(4,k)-Ydoti(k,i))))/rhoi(k,i);
        philower = wrapToPi(-pi/2+thetai(k,i));
        phiupper = wrapToPi(pi/2+thetai(k,i));
        if philower < phiupper
            if ~(phii(k,i) >= philower && phii(k,i) <= phiupper)
                satID(k,i) = NaN;
                rhoi(k,i) = NaN;
                rhoidot(k,i) = NaN;
                phii(k,i) = NaN;
            end
        else
            if ~(phii(k,i) >= philower || phii(k,i) <= phiupper)
                satID(k,i) = NaN;
                rhoi(k,i) = NaN;
                rhoidot(k,i) = NaN;
                phii(k,i) = NaN;
            end
        end
    end
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
plot(tspan,rhoi)
title('Linearized ground stations')
xlabel('Time (s)')
ylabel('\rho^i (m)')
subplot(4,1,2)
plot(tspan,rhoidot)
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