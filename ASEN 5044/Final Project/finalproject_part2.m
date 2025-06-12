clear
clc

% ASEN 5044
% Final Project Part 2
% Questions 5 and 6

%% Constants/provided data/initial conditions
rE = 6378; % km
omegaE = 2*pi/86400; % rad/s
alt = 300; % km
r0 = rE + alt;
mu = 398600; % km^3/s^2
stations = 12;
rng(100);

data = load('orbitdeterm_finalproj_KFdata.mat');
ydata = data.ydata;
Qtrue = data.Qtrue;
Rtrue = data.Rtrue;

% initial conditions
perturb_x0 = [0;0.075;0;-0.021];
x0 = [r0;0;0;r0*sqrt(mu/r0^3)];

dt = 10; % s
tp = 2*pi*sqrt(r0^3/mu);
t0 = 0;
tf = 1400*dt;
tspan = t0:dt:tf;

Nruns = 1;
alpha = 0.05;
NEES = zeros(Nruns,length(tspan));
NIS = zeros(Nruns,length(tspan));

xKF = NaN(Nruns,4,length(tspan));
err = NaN(Nruns,4,length(tspan));
sigma = NaN(Nruns,4,length(tspan));

P0_LKF = 1e-1*diag([.5 .0001 .5 .0001]);
P0_EKF = diag([.5 .1 .5 .1]);

gamma = [0 0; 1 0; 0 0; 0 1];
OmegaKF = dt*gamma;

Q_LKF = 1e-3*[1 -1e-8;-1e-8 1];
Q_EKF = 1e-10*[1 -1e-3;-1e-3 1];

xnom = tmtxdata(x0+perturb_x0,tspan,0,mu,Qtrue);

for n=1:Nruns
    
    %% Simulated ground truth trajectories and measurements using ode45

    KFperturb_x0 = mvnrnd(perturb_x0,P0_LKF)';
    KFx0 = x0+KFperturb_x0;
    

    xsim = tmtxdata(KFx0,tspan,1,mu,Qtrue);
    [xs,ys] = tmtydata(xsim,tspan,1,stations,rE,omegaE,Rtrue);
    [xs_nom,ys_nom] = tmtydata(xnom,tspan,0,stations,rE,omegaE,Rtrue);


    %[x0_LLS,P0_LLS] = warmstartLLS(xnom(:,1:2),xs_nom(:,:,1:2),ys(:,:,1:2),Rtrue);
    
    %% Kalman Filter
    % comment/uncomment to select LKF/EKF/UKF
    % [xKF(n,:,:),err(n,:,:),sigma(n,:,:),NEES(n,:),NIS(n,:),r1x,r2x,r1y,r2y] = LKF(tspan,dt,perturb_x0,P0_LKF,xsim,xnom,xs_nom,ys,OmegaKF,Q_LKF,Rtrue,mu,alpha,Nruns);
    [xKF(n,:,:),err(n,:,:),sigma(n,:,:),NEES(n,:),NIS(n,:),r1x,r2x,r1y,r2y] = EKF(tspan,dt,x0+perturb_x0,P0_EKF,xsim,xs_nom,ys,OmegaKF,Q_EKF,Rtrue,mu,alpha,Nruns);
    %[xKF(n,:,:),err(n,:,:),sigma(n,:,:),NEES(n,:),NIS(n,:),r1x,r2x,r1y,r2y] = UKF(tspan,dt,x0+perturb_x0,P0_EKF,xsim,xs_nom,ys,OmegaKF,Q_EKF,Rtrue,mu,alpha,Nruns);
end

%% NEES/NIS

epsNEES = mean(NEES,1);
epsNIS = mean(NIS,1);


%% Plot ground truth data
figure()
subplot(4,1,1)
plot(tspan,xnom(1,:),tspan,xsim(1,:))
ylabel("X [m]")
legend("Nominal Trajectory", "Simulated Noisy Trajectory")
subplot(4,1,2)
plot(tspan,xnom(2,:),tspan,xsim(2,:))
ylabel("Xdot [m/s]")
subplot(4,1,3)
plot(tspan,xnom(3,:),tspan,xsim(3,:))
ylabel("Y [m]")
subplot(4,1,4)
plot(tspan,xnom(4,:),tspan,xsim(4,:))
xlabel("Time [sec]")
ylabel("Ydot [m/s]")
sgtitle("Noisy simulated ground truth states")

%%

figure()
subplot(4,1,1)
for i=1:stations
    plot(tspan,squeeze(ys(1,i,:)),'+')
    hold on
end
hold off
ylabel('\rho^i (m)')
subplot(4,1,2)
for i=1:stations
    plot(tspan,squeeze(ys(2,i,:)),'+')
    hold on
end
hold off
ylabel('\rho^idot (m/s)')
subplot(4,1,3)
for i=1:stations
    plot(tspan,squeeze(ys(3,i,:)),'+')
    hold on
end
ylabel('\phi^i (radians)')
station_id = NaN(stations, length(tspan));
for i = 1:length(tspan)
    for j = 1:stations
        if ~isnan(ys(1,j,i))
            station_id(j,i) = j;
        end
    end
end
subplot(4,1,4)
plot(tspan, station_id,'+')
sgtitle("Noisy Simulated Measurement Data")
xlabel('Time (s)')
ylabel("Station Access")


%%
figure()
subplot(4,1,1)
plot(tspan,squeeze(err(n,1,:)));
xlabel('time (s)')
ylabel('x_1 (m) error')
subplot(4,1,2)
plot(tspan,squeeze(err(n,2,:)));
xlabel('time (s)')
ylabel('x_2 (m/s) error')
subplot(4,1,3)
plot(tspan,squeeze(err(n,3,:)));
xlabel('time (s)')
ylabel('x_3 error (m)')
subplot(4,1,4)
plot(tspan,squeeze(err(n,4,:)));
xlabel('time (s)')
ylabel('x_4 error (m/s)')
sgtitle('EKF errors')

%%

figure()
subplot(2,2,1)
plot(tspan,squeeze(err(n,1,:)),tspan,2*squeeze(sigma(n,1,:)),'r',tspan,-2*squeeze(sigma(n,1,:)),'r');
xlabel('time (s)')
ylabel('x_1 (m) error')
ylim([-1, 1])
subplot(2,2,2)
plot(tspan,squeeze(err(n,2,:)),tspan,2*squeeze(sigma(n,2,:)),'r',tspan,-2*squeeze(sigma(n,2,:)),'r');
xlabel('time (s)')
ylabel('x_2 (m/s) error')
ylim([-0.01, 0.01])
subplot(2,2,3)
plot(tspan,squeeze(err(n,3,:)),tspan,2*squeeze(sigma(n,3,:)),'r',tspan,-2*squeeze(sigma(n,3,:)),'r');
xlabel('time (s)')
ylabel('x_3 error (m)')
ylim([-1, 1])
subplot(2,2,4)
plot(tspan,squeeze(err(n,4,:)),tspan,2*squeeze(sigma(n,4,:)),'r',tspan,-2*squeeze(sigma(n,4,:)),'r');
xlabel('time (s)')
ylabel('x_4 error (m/s)')
ylim([-0.01, 0.01])
legend("Error","2\sigma bound")
sgtitle('EKF errors')

%%
figure()
subplot(4,1,1)
plot(tspan,squeeze(xKF(n,1,:)));
xlabel('time (s)')
ylabel('x_1 (m)')
subplot(4,1,2)
plot(tspan,squeeze(xKF(n,2,:)));
xlabel('time (s)')
ylabel('x_2 (m/s)')
subplot(4,1,3)
plot(tspan,squeeze(xKF(n,3,:)));
xlabel('time (s)')
ylabel('x_3 (m)')
subplot(4,1,4)
plot(tspan,squeeze(xKF(n,4,:)));
xlabel('time (s)')
ylabel('x_4 (m/s)')
%sgtitle('EKF States')
%%
figure()
subplot(2,1,1)
for i=1:length(Nruns)
    plot(tspan,epsNEES(i,:),'+',tspan,r1x,tspan,r2x)
    hold on
end
ylim([-10 20])
title('NEES')
ylabel('NEES statistic, $\bar{\epsilon}_x$','interpreter', 'latex', 'FontSize',14)
xlabel('time step, k','FontSize',14)
title('NEES Estimation Results','FontSize',14)
ylim([2, 6])
legend('NEES @ time k', 'r_1 bound', 'r_2 bound') 
subplot(2,1,2)
for i=1:length(Nruns)
    plot(tspan,epsNIS(i,:),'+',tspan,r1y,tspan,r2y)
    hold on
end
title('NIS')
%sgtitle('EKF NEES & NIS')
ylim([-2 10])
ylabel('NIS statistic, $\bar{\epsilon}_y$','interpreter', 'latex','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')

%% Calculate solution using "real" data
ylog = NaN(3,12,length(tspan));
for i=1:length(tspan)
    y = cell2mat(ydata(i));
    if ~isnan(y)
        for j=1:length(y(1,:))
            stationID = y(4,j);
            ylog(:,stationID,i) = y(1:3,j);
        end
    end

end

% Plot ground stations
figure()
subplot(3,1,1)
plot(tspan,squeeze(ylog(1,:,:)), 'Marker','x')
xlabel('Time (s)')
ylabel('\rho^i (m)')
subplot(3,1,2)
plot(tspan,squeeze(ylog(2,:,:)), 'Marker','o')
xlabel('Time (s)')
ylabel('\rho^idot (m/s)')
subplot(3,1,3)
plot(tspan,squeeze(ylog(3,:,:)),'+')
xlabel('Time (s)')
ylabel('\phi^i (radians)')
sgtitle("Given Observation Data")

% Comment to select LKF/EKF for filtering observation data
[xKFobs,err_obs,sigma_obs,NEESobs,NISobs,r1x_obs,r2x_obs,r1y_obs,r2y_obs] = LKF(tspan,dt,perturb_x0,P0_LKF,xsim,xnom,xs_nom,ylog,OmegaKF,Q_LKF,Rtrue,mu,alpha,Nruns);
% [xKFobs,err_obs,sigma_obs,NEESobs,NISobs,r1x_obs,r2x_obs,r1y_obs,r2y_obs] = EKF(tspan,dt,x0,P0_EKF,xsim,xs_nom,ylog,OmegaKF,Q_EKF,Rtrue,mu,alpha,Nruns);

% Plot observation data
figure()
subplot(4,1,1)
plot(tspan,err_obs(1,:));
subplot(4,1,2)
plot(tspan,err_obs(2,:));
subplot(4,1,3)
plot(tspan,err_obs(3,:));
subplot(4,1,4)
plot(tspan,err_obs(4,:));


figure()
subplot(4,1,1)
plot(tspan,err_obs(1,:),tspan,err_obs(1,:)+2*sigma_obs(1,:),'r',tspan,err_obs(1,:)-2*sigma_obs(1,:),'r');
xlabel('time (s)')
ylabel('x_1 (m) error')
% ylim([-0.8, 0.8])
subplot(4,1,2)
plot(tspan,err_obs(2,:),tspan,err_obs(2,:)+2*sigma_obs(2,:),'r',tspan,err_obs(2,:)-2*sigma_obs(2,:),'r');
xlabel('time (s)')
ylabel('x_2 (m/s) error')
% ylim([-0.01 0.01])
subplot(4,1,3)
plot(tspan,err_obs(3,:),tspan,err_obs(3,:)+2*sigma_obs(3,:),'r',tspan,err_obs(3,:)-2*sigma_obs(3,:),'r');
xlabel('time (s)')
ylabel('x_3 error (m)')
% ylim([-0.8, 0.8])
subplot(4,1,4)
plot(tspan,err_obs(4,:),tspan,err_obs(4,:)+2*sigma_obs(4,:),'r',tspan,err_obs(4,:)-2*sigma_obs(4,:),'r');
xlabel('time (s)')
ylabel('x_4 error (m/s)')
% ylim([-0.01 0.01])
legend("Error", "2\sigma bound")
sgtitle('LKF errors')

figure()
subplot(4,1,1)
plot(tspan,xKFobs(1,:));
xlabel('time (s)')
ylabel('x_1 (m)')
subplot(4,1,2)
plot(tspan,xKFobs(2,:));
xlabel('time (s)')
ylabel('x_2 (m/s)')
subplot(4,1,3)
plot(tspan,xKFobs(3,:));
xlabel('time (s)')
ylabel('x_3 (m)')
subplot(4,1,4)
plot(tspan,xKFobs(4,:));
xlabel('time (s)')
ylabel('x_4 (m/s)')
sgtitle('EKF States')


%% Support functions

% nonlinear system
function xdot = eom(x,mu)
    r = sqrt(x(1)^2+x(3)^2);
    xdot = [x(2); -mu*x(1)/r^3; x(4); -mu*x(3)/r^3];
end

function xdot = eom_noise(x,mu,noise)
    r = sqrt(x(1)^2+x(3)^2);
    xdot = [x(2); -mu*x(1)/r^3 + noise(1); x(4); -mu*x(3)/r^3 + noise(2)];
end

% measurements

function xdata = tmtxdata(x0,tspan,noiseYN,mu,Qtrue)
    tol = 1e-12;
    options = odeset('Stats','off','RelTol',tol,'AbsTol',tol);
    xdata = NaN(length(x0),length(tspan));
    xdata(:,1) = x0;
    dt = tspan(2)-tspan(1);
    for k=2:length(tspan)
        if noiseYN == 1
            noise = mvnrnd([0 0],Qtrue)';
        else
            noise = [0;0];
        end
        [tout, xout] = ode45(@(tout,xout) eom_noise(xout,mu,noise),[0 dt],xdata(:,k-1),options);
        xdata(:,k) = xout(end,:)';
    end
end

function [xs,ys] = tmtydata(xdata,tspan,noiseYN,stations,rE,omegaE,Rtrue)
    % measurements
    xs = NaN(4,stations,length(tspan));
    ys = NaN(3,stations,length(tspan));

    % simulated measurements
    for k=1:length(tspan)
        for i=1:stations
            thetai0 = (i-1)*pi/6;
            xs(1,i,k) = rE*cos(omegaE*tspan(k)+thetai0);
            xs(3,i,k) = rE*sin(omegaE*tspan(k)+thetai0);
            xs(2,i,k) = -omegaE*xs(3,i,k);
            xs(4,i,k) = omegaE*xs(1,i,k);
            if k>1
                thetai = atan2(xs(3,i,k),xs(1,i,k));
            else
                thetai = thetai0;
            end
    
            ys(:,i,k) = hmatrix(xdata(:,k),xs(:,i,k));
            %ys(:,i,k) = ys(:,i,k) + noise(:,k);
    
            philower = wrapToPi(-pi/2+thetai);
            phiupper = wrapToPi(pi/2+thetai);

            if noiseYN==1
                noise = mvnrnd([0 0 0],Rtrue)';
            else
                noise = [0;0;0];
            end

            if philower < phiupper
                if ~(ys(3,i,k) >= philower && ys(3,i,k) <= phiupper)
                    ys(:,i,k) = NaN;
                else
                    ys(:,i,k) = ys(:,i,k) + noise;
                end
            else
                if ~(ys(3,i,k) >= philower || ys(3,i,k) <= phiupper)
                    ys(:,i,k) = NaN;
                else
                    ys(:,i,k) = ys(:,i,k) + noise;
                end
            end
        end
    end
end

% initial conditions

function [mu0, p0] = warmstartLLS(xdata,xs,ys,Rtrue)
    time = length(ys(1,1,:));
    %R = repmat({Rtrue},time,1);
    %R = blkdiag(R{:});
    %H = NaN(3*time,4);
    R = Rtrue;
    H = [];
    bigystack = [];
    smallH = [];
    dy = [];
    for k=1:time
        StationID = find(~isnan(squeeze(ys(1,:,k))));
        ystack = NaN(3*length(StationID),1);
            for i=1:length(StationID)
                ystack(1+3*(i-1):3+3*(i-1),1) = ys(:,StationID(i),k);
                smallH = [smallH; Cmatrix(xdata(:,k),xs(:,StationID(i),k))];
                R = blkdiag(R,Rtrue);
            end
        bigystack = [bigystack; ystack];
        %H(1+3*(k-1):3+3*(k-1),:) = smallH;
        H = [H; smallH];
        dy = [dy; bigystack - smallH*xdata(:,k)];
        %R = blkdiag(R,Rtrue);
        %dy(1+3*(k-1):3+3*(k-1)) = ystack(1+3*(k-1):3+3*(k-1),1) - smallH*xdata(:,k);
    end
    mu0 = (H'*R^(-1)*H)^(-1)*H'*R^(-1)*dy;
    p0 = inv(H'*inv(R)*H);
end

% matrices/jacobians

function A = Amatrix(x,mu)
    A  = [0 1 0 0; ...
        -mu*(x(1)^2+x(3)^2)^(-3/2) + 3*mu*x(1)^2*(x(1)^2+x(3)^2)^(-5/2),...
        0,...
        3*mu*x(1)*x(3)*(x(1)^2+x(3)^2)^(-5/2),...
        0; ...
        0 0 0 1;...
        3*mu*x(1)*x(3)*(x(1)^2+x(3)^2)^(-5/2),...
        0,...
        -mu*(x(1)^2+x(3)^2)^(-3/2) + 3*mu*x(3)^2*(x(1)^2+x(3)^2)^(-5/2),...
        0];
end

function C = Cmatrix(xnom,xi)
    Rs = sqrt((xnom(1)-xi(1))^2 + (xnom(3)-xi(3))^2);
    rho_n = (xnom(1)-xi(1))*(xnom(2)-xi(2)) + (xnom(3)-xi(3))*(xnom(4)-xi(4));
    C = [(xnom(1)-xi(1))/Rs, 0, (xnom(3)-xi(3))/Rs, 0; ...
        (xnom(2)-xi(2))/Rs - (xnom(1)-xi(1))*rho_n/Rs^3, (xnom(1)-xi(1))/Rs, (xnom(4)-xi(4))/Rs - (xnom(3)-xi(3))*rho_n/Rs^3, (xnom(3)-xi(3))/Rs; ...
        -(xnom(3)-xi(3))/Rs^2, 0, (xnom(1)-xi(1))/Rs^2, 0];
end

function y = hmatrix(x,xs)
    % rho
    y1 = sqrt((x(1)-xs(1))^2 + (x(3)-xs(3))^2);
    % rho dot
    y2 = ((x(1)-xs(1))*(x(2)-xs(2)) + (x(3)-xs(3))*(x(4)-xs(4)))/y1;
    % phi
    y3 = atan2((x(3)-xs(3)),x(1)-xs(1));
    y = [y1;y2;y3];
end

% Kalman filters

function [xKF,err,sigma,NEES,NIS,r1x,r2x,r1y,r2y] = LKF(tspan,dt,dx0,P0,xsim,xnom,xs_nom,ys,OmegaKF,Q_KF,Rtrue,mu,alpha,N)

    dxm = NaN(4,length(tspan));
    dxp = NaN(4,length(tspan));
    xKF = NaN(4,length(tspan));
    Pm = NaN(4,4,length(tspan));
    Pp = NaN(4,4,length(tspan));
    err = NaN(4,length(tspan));
    sigma = NaN(4,length(tspan));
    NEES = zeros(1,length(tspan));
    NIS = zeros(1,length(tspan));
    
    dxp(:,1) = dx0;
    Pp(:,:,1) = P0;
    sigma(:,1) = sqrt(diag(P0));

    for k=1:length(tspan)-1
        % A&F matrices
        A = Amatrix(xnom(:,k),mu);
        F = eye(4) + dt*A;
        % Predictor
        dxm(:,k+1) = F*dxp(:,k);
        Pm(:,:,k+1) = F*Pp(:,:,k)*F' + OmegaKF*Q_KF*OmegaKF';
        % Corrector
        if nnz(ys(:,:,k+1))>0
            numStations = sum(~isnan(squeeze(ys(1,:,k+1))));
            StationID = find(~isnan(squeeze(ys(1,:,k+1))));
            ystack = NaN(3*numStations,1);
            ynomstack = NaN(3*numStations,1);
            R = repmat({Rtrue},numStations,1);
            R = blkdiag(R{:});
            H = NaN(length(ystack(:,1)),4);
            for i=1:numStations
                ystack(1+3*(i-1):3+3*(i-1),1) = ys(:,StationID(i),k+1);
                ynomstack(1+3*(i-1):3+3*(i-1),1) = hmatrix(xnom(:,k+1),xs_nom(:,StationID(i),k+1));
                H(1+3*(i-1):3+3*(i-1),:) = Cmatrix(xnom(:,k+1),xs_nom(:,StationID(i),k+1));
            end
            %dy = ystack - H*dxm(:,k+1);
            dy = ystack - ynomstack;
            K = Pm(:,:,k+1)*H'/(H*Pm(:,:,k+1)*H' + R);
            dxp(:,k+1) = dxm(:,k+1) + K*(dy - H*dxm(:,k+1));
            Pp(:,:,k+1) = (eye(4) - K*H)*Pm(:,:,k+1);
            S = H*Pm(:,:,k+1)*H' + R;
            S = 0.5*(S+S');
        else
            dxp(:,k+1) = dxm(:,k+1);
            Pp(:,:,k+1) = Pm(:,:,k+1);
            S = Pm(:,:,k+1) + Rtrue;
            S = 0.5*(S+S');
            dy = [0;0;0];
        end
        xKF(:,k+1) = xnom(:,k+1) + dxp(:,k+1);
        %err(:,k+1) = dxp(:,k+1)-dxm(:,k+1);
        err(:,k+1) = xsim(:,k+1) - xKF(:,k+1);
        sigma(:,k+1) = sqrt(diag(Pp(:,:,k+1)));
        NEES(k+1) = (xsim(:,k+1)-xKF(:,k+1))'/Pp(:,:,k+1)*(xsim(:,k+1)-xKF(:,k+1));
        NIS(k+1) = dy'/S*dy;
        r1x = chi2inv(alpha/2,N*length(err(:,k+1)))./N.*ones(1,length(tspan));
        r2x = chi2inv(1-alpha/2,N*length(err(:,k+1)))./N.*ones(1,length(tspan));
        r1y = chi2inv(alpha/2,N*length(dy))./N.*ones(1,length(tspan));
        r2y = chi2inv(1-alpha/2,N*length(dy))./N.*ones(1,length(tspan));
    end
end

function [xKF,err,sigma,NEES,NIS,r1x,r2x,r1y,r2y] = EKF(tspan,dt,x0,P0,xsim,xs,ys,OmegaKF,Q_KF,Rtrue,mu,alpha,N)

    xm = NaN(4,length(tspan));
    xp = NaN(4,length(tspan));
    xKF = NaN(4,length(tspan));
    Pm = NaN(4,4,length(tspan));
    Pp = NaN(4,4,length(tspan));
    err = NaN(4,length(tspan));
    sigma = NaN(4,length(tspan));
    NEES = zeros(1,length(tspan));
    NIS = zeros(1,length(tspan));
    
    xp(:,1) = x0;
    Pp(:,:,1) = P0;
    sigma(:,1) = sqrt(diag(P0));
    
    tol = 1e-12;
    options = odeset('Stats','off','RelTol',tol,'AbsTol',tol);

    for k=1:length(tspan)-1
        % AC matrices
        A = Amatrix(xp(:,k),mu);
        F = eye(4) + dt*A;
        % Predictor
        [tout, xout] = ode45(@(tout,xout) eom(xout,mu),[0 dt],xp(:,k),options);
        xm(:,k+1) = xout(end,:)';
        Pm(:,:,k+1) = F*Pp(:,:,k)*F' + OmegaKF*Q_KF*OmegaKF';
        % Corrector
        if nnz(ys(:,:,k+1))>0
            numStations = sum(~isnan(squeeze(ys(1,:,k+1))));
            StationID = find(~isnan(squeeze(ys(1,:,k+1))));
            ystack = NaN(3*numStations,1);
            ymstack = NaN(3*numStations,1);
            R = repmat({Rtrue},numStations,1);
            R = blkdiag(R{:});
            H = NaN(length(ystack(:,1)),4);
            for i=1:numStations
                ystack(1+3*(i-1):3+3*(i-1),1) = ys(:,StationID(i),k+1);
                ymstack(1+3*(i-1):3+3*(i-1),1) = hmatrix(xm(:,k+1),xs(:,StationID(i),k+1));
                H(1+3*(i-1):3+3*(i-1),:) = Cmatrix(xm(:,k+1),xs(:,StationID(i),k+1));
            end
            dy = ystack - ymstack;
            K = Pm(:,:,k+1)*H'/(H*Pm(:,:,k+1)*H' + R);
            xp(:,k+1) = xm(:,k+1) + K*dy;
            Pp(:,:,k+1) = (eye(4) - K*H)*Pm(:,:,k+1);
            S = H*Pm(:,:,k+1)*H' + R;
            S = 0.5*(S+S');
        else
            xp(:,k+1) = xm(:,k+1);
            Pp(:,:,k+1) = Pm(:,:,k+1);
            S = Pm(:,:,k+1) + Rtrue;
            S = 0.5*(S+S');
            dy = [0;0;0];
        end
        xKF(:,k+1) = xp(:,k+1);
        err(:,k+1) = xp(:,k+1) - xm(:,k+1);
        sigma(:,k+1) = sqrt(diag(Pp(:,:,k+1)));
        NEES(k+1) = (xsim(:,k+1)-xp(:,k+1))'*inv(Pp(:,:,k+1))*(xsim(:,k+1)-xp(:,k+1));
        NIS(k+1) = dy'*inv(S)*dy;
        r1x = chi2inv(alpha/2,N*length(err(:,k+1)))./N.*ones(1,length(tspan));
        r2x = chi2inv(1-alpha/2,N*length(err(:,k+1)))./N.*ones(1,length(tspan));
        r1y = chi2inv(alpha/2,N*length(dy))./N.*ones(1,length(tspan));
        r2y = chi2inv(1-alpha/2,N*length(dy))./N.*ones(1,length(tspan));
    end
end

function [xKF,err,sigma,NEES,NIS,r1x,r2x,r1y,r2y] = UKF(tspan,dt,x0,P0,xsim,xs,ys,OmegaKF,Q_KF,Rtrue,mu,alphaN,N)

    tol = 1e-12;
    options = odeset('Stats','off','RelTol',tol,'AbsTol',tol);

    xp = NaN(4,length(tspan));
    xm = NaN(4,length(tspan));
    Pp = NaN(4,4,length(tspan));
    Pm = NaN(4,4,length(tspan));
    Pyy = NaN(3,3,length(tspan));
    Pxy = zeros(4,3,length(tspan));
    xKF = NaN(4,length(tspan));
    err = NaN(4,length(tspan));
    sigma = NaN(4,length(tspan));
    NEES = zeros(1,length(tspan));
    NIS = zeros(1,length(tspan));

    xp(:,1) = x0;
    Pp(:,:,1) = P0;
    sigma(:,1) = sqrt(diag(P0));

    n = length(x0);
    kappa = 0;
    beta = 2;
    alpha = 1e-3;
    lambda = alpha^2*(n+kappa)-n;

    wm = NaN(2*n+1,1);
    wc = NaN(2*n+1,1);
    wm(1) = lambda/(n+lambda);
    wc(1) = lambda/(n+lambda) + 1 - alpha^2 + beta;
    for i=1:2*n
        wm(i+1) = 1/(2*(n+lambda));
        wc(i+1) = wm(i+1);
    end

    for k = 1:length(tspan)-1
        chi = NaN(n,2*n+1);
        chi(:,1) = xp(:,k);
        S = chol(Pp(:,:,k));
        xm(:,k+1) = zeros(4,1);
        Pm(:,:,k+1) = OmegaKF*Q_KF*OmegaKF';
        for i=1:2*n
            if i<=n
                chi(:,i+1) = xp(:,k) + sqrt(n+lambda)*S(i,:)';
            else
                chi(:,i+1) = xp(:,k) - sqrt(n+lambda)*S(i-n,:)';
            end
        end
        for i=1:2*n+1
            [tout xout] = ode45(@(tout,xout) eom(xout,mu), [0 dt], chi(:,i), options);
            chibar = xout(end,:)';
            xm(:,k+1) = xm(:,k+1) + wm(i)*chibar;
            Pm(:,:,k+1) = Pm(:,:,k+1) + wc(i)*(chibar-xm(:,k+1))*(chibar-xm(:,k+1))';
        end

        if nnz(ys(:,:,k+1))>0
            numStations = sum(~isnan(squeeze(ys(1,:,k+1))));
            StationID = find(~isnan(squeeze(ys(1,:,k+1))));
            ystack = NaN(3*numStations,1);
            ymstack = NaN(3*numStations,1);
            R = repmat({Rtrue},numStations,1);
            R = blkdiag(R{:});
            H = NaN(length(ystack(:,1)),4);
            for i=1:numStations
                ystack(1+3*(i-1):3+3*(i-1),1) = ys(:,StationID(i),k+1);
                %ymstack(1+3*(i-1):3+3*(i-1),1) = hmatrix(xm(:,k+1),xs(:,StationID(i),k+1));
                H(1+3*(i-1):3+3*(i-1),:) = Cmatrix(xm(:,k+1),xs(:,StationID(i),k+1));
            end
            S = chol(Pm(:,:,k+1));
            chi(:,1) = xm(:,k+1);
            ym = zeros(3*numStations,1);
            Pyy(:,:,k+1) = Rtrue;
            for i=1:2*n
                if i<=n
                    chi(:,i+1) = xm(:,k+1) + sqrt(n+lambda)*S(i,:)';
                else
                    chi(:,i+1) = xm(:,k+1) - sqrt(n+lambda)*S(i-n,:)';
                end
            end
            for i=1:2*n+1
                gamma = zeros(3*numStations,1);
                for s=1:numStations
                    gamma(1+3*(s-1):3+3*(s-1),1) = hmatrix(chi(:,i),xs(:,StationID(s),k+1));
                end
                ym = ym + wm(i)*gamma;
                Pyy(:,:,k+1) = Pyy(:,:,k+1) + wc(i)*(gamma-ym)*(gamma-ym)';
                Pxy(:,:,k+1) = Pxy(:,:,k+1) + wc(i)*(chi(:,k+1)-xm(:,k+1))*(gamma-ym)';
            end
            K = Pxy(:,:,k+1)/Pyy(:,:,k+1);
            dy = ys(:,k+1) - ym;
            xp(:,k+1) = xm(:,k+1) + K*dy;
            Pp(:,:,k+1) = Pm(:,:,k+1) - K*Pyy(:,:,k+1)*K';
        else
            xp(:,k+1) = xm(:,k+1);
            Pp(:,:,k+1) = Pm(:,:,k+1);
        end
        
        xKF(:,k+1) = xp(:,k+1);
        err(:,k+1) = xsim(:,k+1) - xKF(:,k+1);
        sigma(:,k+1) = sqrt(diag(Pp(:,:,k+1)));
        % NEES(k+1) = (err(:,k+1))'*inv(Pp(:,:,k+1))*(err(:,k+1));
        % NIS(k+1) = dy'*inv(S)*dy;
        % r1x = chi2inv(alphaN/2,N*length(err(:,k+1)))./N.*ones(1,length(tspan));
        % r2x = chi2inv(1-alphaN/2,N*length(err(:,k+1)))./N.*ones(1,length(tspan));
        % r1y = chi2inv(alphaN/2,N*length(dy))./N.*ones(1,length(tspan));
        % r2y = chi2inv(1-alphaN/2,N*length(dy))./N.*ones(1,length(tspan));
    end
end