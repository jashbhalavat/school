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
% x0 = [r0; 0; 0; r0*sqrt(mu/r0^3)];
x0 = [r0; 0; 0; r0*sqrt(mu/r0^3)] + perturb_x0;

dt = 10; % s
t0 = 0;
tp = 2*pi*sqrt(r0^3/mu);
tf = 1400*dt;

tspan = t0:dt:tf;

% ode45 solution
tol = 1e-12;
options = odeset('RelTol',tol);
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
   satID(1,i) = i; 
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

   % sc_from_gs = [xout(1,1) - Xi(1,i), xout(1,3) - Yi(1,i)];
   % angle = ang_diff(sc_from_gs, [Xi(1,i), Yi(1,i)]);
   % if angle > abs(pi/2)
   %     satID(1,i) = NaN;
   %     rhoi(1,i) = NaN;
   %     rhoidot(1,i) = NaN;
   %     phii(1,i) = NaN;
   % end
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
        xi = [Xi(k,i), Xdoti(k,i), Yi(k,i), Ydoti(k,i)]';
        xd = xout(k,:) - xi';
        thetai(k,i) = atan2(Yi(k,i),Xi(k,i));
        phii(k,i) = atan2(xd(3),xd(1));
        rhoi(k,i) = sqrt(xd(1)^2 + xd(3)^2);
        rhoidot(k,i) = (xd(1)*xd(2) + xd(3)*xd(4))/rhoi(k,i);
        philower = wrapToPi(-pi/2+thetai(k,i));
        phiupper = wrapToPi(pi/2+thetai(k,i));

        % sc_from_gs = [xout(k,1) - Xi(k,i), xout(k,3) - Yi(k,i)];
        % angle = ang_diff(sc_from_gs, [Xi(k,i), Yi(k,i)]);
        % if angle > abs(pi/2)
        %    satID(k,i) = NaN;
        %    rhoi(k,i) = NaN;
        %    rhoidot(k,i) = NaN;
        %    phii(k,i) = NaN;
        % end

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
plot(tout,rhoi, 'Marker','x')
title('Tracking stations with ode45')
xlabel('Time (s)')
ylabel('\rho^i (m)')
subplot(4,1,2)
plot(tout,rhoidot, 'Marker','o')
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

function angle = ang_diff(V1, V2)
    angle = acos(dot(V1, V2) / (norm(V1) * norm(V2)));
end

