clear; clc; close all;

% ASEN 5044 - Final Project - Linearized KF

%% Constants

rng(100);
r0 = 6678; % km
rE = 6378; % km
mu = 398600; % km^3/s^2
k = sqrt(mu/r0^3); % mean motion
dt = 10; % sec
omegaE = 2*pi/86400; % rad/s
period = 2*pi*sqrt(r0^3/mu); % s
stations = 12;
alt = 300; % km
r0 = rE + alt;
num_sims = 5;

% System variables
n = 4;
p = 3;
I_n = eye(n);
I_p = eye(p);

% Initial conditions
perturb_x0 = [0;0.075;0;-0.021];
x0 = [r0; 0; 0; r0*sqrt(mu/r0^3)] + perturb_x0;

dt = 10; % s
t0 = 0; % s
tp = 2*pi*sqrt(r0^3/mu); % s
tf = 1400*dt; % s
tspan = t0:dt:tf; % s

%% Parse through given data

data = load("orbitdeterm_finalproj_KFdata.mat");
Qtrue = data.Qtrue;
Rtrue = data.Rtrue;
measLabels = data.measLabels;
tvec = data.tvec;
raw_y_data = data.ydata;

rho_i_true = NaN(length(tvec), stations); % km
rho_dot_i_true = NaN(length(tvec), stations); % km/s
phi_i_true = NaN(length(tvec), stations); % rad
station_id_true = NaN(length(tvec), stations);

for i = 1:length(raw_y_data)
    cell_content = raw_y_data{i};
    cell_cols = size(cell_content, 2);
    for j = 1:cell_cols
        if ~isnan(cell_content(4, j))
            rho_i_true(i, cell_content(4, j)) = cell_content(1, j);
            rho_dot_i_true(i, cell_content(4, j)) = cell_content(2, j);
            phi_i_true(i, cell_content(4, j)) = cell_content(3, j);
            station_id_true(i, cell_content(4, j)) = cell_content(4, j);
        end
    end
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
xi = NaN(4, stations, length(tspan));

for k = 1:length(tspan)
    for i=1:stations
       satID_nom(k,i) = i;
       thetai(1,i) = (i-1)*pi/6;
       Xi(k,i) = rE*cos(omegaE*tspan(k) + thetai(1,i));
       Yi(k,i) = rE*sin(omegaE*tspan(k) + thetai(1,i));
       Xdoti(k,i) = -omegaE*Yi(k,i);
       Ydoti(k,i) = omegaE*Xi(k,i);
       xi(:,i,k) = [Xi(k,i); Xdoti(k,i); Yi(k,i); Ydoti(k,i)];
       phii_nom(k,i) = atan2(xout(k,3)-Yi(k,i), xout(k,1)-Xi(k,i));
       rhoi_nom(k,i) = sqrt((xout(k,1)-Xi(k,i))^2 + (xout(k,3)-Yi(k,i))^2);
       rhoidot_nom(k,i) = ((xout(k,1)-Xi(k,i))*(xout(k,2)-Xdoti(k,i)) + ((xout(k,3)-Yi(k,i))*(xout(k,4)-Ydoti(k,i))))/rhoi_nom(k,i);

       if k > 1
            thetai(k,i) = atan2(Yi(k,i),Xi(k,i));
       end
       
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

A(:,:,1) = A_func(xnom(:,1), mu);
F(:,:,1) = eye(4) + A(:,:,1)*dt;

for i=2:length(tspan)
    % Nominal state
    x = r0*cos(2*pi/tp*tspan(i));
    y = r0*sin(2*pi/tp*tspan(i));
    vx = -2*pi/tp*y;
    vy = 2*pi/tp*x;
    xnom(:,i) = [x; vx; y; vy];

    % perturbations linearized about nominal point
    A(:,:,i) = A_func(xnom(:,i), mu);
    F(:,:,i) = eye(4) + dt * A(:,:,i);
    dx(:,i) = F(:,:,i-1)*dx(:,i-1);
    xtot(:,i) = xnom(:,i) + dx(:,i);
end

phii = NaN(length(tspan),stations);
rhoi = NaN(length(tspan),stations);
rhoidot = NaN(length(tspan),stations);
satID = NaN(length(tspan),stations);
ynom = NaN(3,stations, length(tspan));
ytot = NaN(3,stations, length(tspan));
dy = NaN(3,stations, length(tspan));

for k = 1:length(tspan)
    for i=1:stations
       satID(k,i) = i;
       
       ynom(:,i,k) = h(xnom(:,k), xi(:,i,k));
       dhdx = dh_dx(ynom(:,i,k), xi(:,i,k), xnom(:,k));
       dy(:,i,k) = dhdx * dx(:,k);
       ytot(:,i,k) = ynom(:,i,k) + dy(:,i,k);
       rhoi(k,i) = ytot(1,i,k);
       rhoidot(k,i) = ytot(2,i,k);
       phii(k,i) = ytot(3,i,k);
       
       if k > 1
            thetai(k,i) = atan2(xi(3,i,k),xi(1,i,k));
       end
        
       if isnan(satID_nom(k, i))
            satID(k,i) = NaN;
            rhoi(k,i) = NaN;
            rhoidot(k,i) = NaN;
            phii(k,i) = NaN;
       end
   end
end

function out = A_func(xnom, mu)
    out = [0 1 0 0; -mu*(xnom(1)^2+xnom(3)^2)^(-3/2) + 3*mu*xnom(1)^2*(xnom(1)^2+xnom(3)^2)^(-5/2) 0 3*mu*xnom(1)*xnom(3)*(xnom(1)^2+xnom(3)^2)^(-5/2) 0; ...
        0 0 0 1; 3*mu*xnom(1)*xnom(3)*(xnom(1)^2+xnom(3)^2)^(-5/2) 0 -mu*(xnom(1)^2+xnom(3)^2)^(-3/2) + 3*mu*xnom(3)^2*(xnom(1)^2+xnom(3)^2)^(-5/2) 0];
end

function out = h(x_star, x_i)
    out(1,1) = sqrt((x_star(1)-x_i(1))^2 + (x_star(3)-x_i(3))^2);
    out(2,1) = ((x_star(1)-x_i(1))*(x_star(2)-x_i(2)) + (x_star(3)-x_i(3))*(x_star(4)-x_i(4)))/out(1);
    out(3,1) = atan2((x_star(3)-x_i(3)),(x_star(1)-x_i(1)));
end

function out = dh_dx(y_nom, X_gs, xnom)
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
    out = [dh1_dx1, 0, dh1_dx3, 0;
        dh2_dx1, dh2_dx2, dh2_dx3, dh2_dx4;
        dh3_dx1, 0, dh3_dx3, 0];
end

%% Monte Carlo Simulation

Qlkf = Qtrue;
Rlkf = blkdiag(Rtrue, Rtrue, Rtrue, Rtrue, Rtrue, Rtrue, Rtrue, Rtrue, Rtrue, Rtrue, Rtrue, Rtrue);
NEESsamps = zeros(num_sims,length(tspan));
NISsamps = zeros(num_sims,length(tspan));

gamma = [0 0; 1 0; 0 0; 0 1];
omega = dt * gamma;

xout = xout';

for ss = 1:num_sims
    % Generate true trajectory and measurements from system
    xk_truehist = zeros(n,length(tspan));
    ykhist = zeros(p,stations,length(tspan));
    
    wk = mvnrnd(zeros(2,1),Qtrue)';
    xk_truehist(:,1) = x0-perturb_x0 + gamma*wk;
    vk = mvnrnd(zeros(p,1),Rtrue)';
    for station = 1:stations
        if ~isnan(satID_nom(1, station))
                ykhist(:,station,1) = h(xk_truehist(:,1), xi(:,station,1)) + vk;
        end
    end

    for jj = 2:length(tspan)
        ti = tspan(jj-1);
        tf = tspan(jj);
        wk = mvnrnd(zeros(2,1),Qtrue)';
        xk = xk_truehist(:,jj-1);
        [tode_out, xode_out] = ode45(@(tout, xk) eom_full_nl(xk,mu,wk),ti:dt/2:tf,xk_truehist(:,jj-1),options);
        xk_truehist(:,jj) = xode_out(end,:)';
        vk = mvnrnd(zeros(p,1),Rtrue)';
        for station = 1:stations
            if ~isnan(satID_nom(jj, station))
                ykhist(:,station,jj) = h(xk_truehist(:,jj), xi(:,station,jj)) + vk;
            end
        end
    end

    % Linearized Kalman Filter
    NEESsshist = zeros(1, length(tspan));
    NISsshist = zeros(1, length(tspan));

    F_lkf = NaN(n,n,length(tspan));
    H_lkf = zeros(p*n,n,length(tspan));
    dy_lkf = zeros(p*n,length(tspan));

    for jj = 1:length(tspan)
        F_lkf(:,:,jj) = I_n + dt*A_func(xout(:,jj), mu);
        for station = 1:stations
            if ~isnan(satID_nom(jj, station))
                H_lkf((3*station-2):(3*station),:,jj) = dh_dx(ykhist(:,station,jj), xi(:,station,jj), xout(:,jj));
                dy_lkf(3*station-2:3*station,jj) = ykhist(:,jj) - h(xout(:,jj), xi(:,station,jj));
            end
        end
    end

    P_plus_0 = eye(n);

    [xlkf_out, Plkf_out, Pyylkf_out, innovlkf_out] = lkf(tspan, xout, F_lkf, perturb_x0, Rlkf, dy_lkf, H_lkf, P_plus_0, Qlkf, omega);

    % Compute and store NEES and NIS statistics:
    for jj = 1:length(tspan)
        invPkp1 = inv(Plkf_out(:,:,jj));
        invPyykp1 = inv(Pyylkf_out(:,:,jj));
        NEESsshist(jj) = ...
            (xk_truehist(:,jj) - xlkf_out(:,jj))'*invPkp1*(xk_truehist(:,jj) - xlkf_out(:,jj));
        NISsshist(jj) = innovlkf_out(:,jj)'*invPyykp1*innovlkf_out(:,jj);
    end
    NEESsamps(ss,:) = NEESsshist;
    NISsamps(ss,:) = NISsshist;

end

%%DO NEES TEST:
epsNEESbar = mean(NEESsamps,1);
alphaNEES = 0.05; %%significance level
Nnx = num_sims*length(F_lkf);
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ num_sims
r2x = chi2inv(1-alphaNEES/2, Nnx )./ num_sims
figure()
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, \bar{\epsilon}_x','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')

%%
%%DO NIS TEST:
epsNISbar = mean(NISsamps,1);
alphaNIS = 0.05; 
Nny = num_sims*size(H_lkf,1);
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ num_sims
r2y = chi2inv(1-alphaNIS/2, Nny )./ num_sims
figure(90)
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, \bar{\epsilon}_y','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')

function xdot = eom_full_nl(x,mu,wk)
    r = sqrt(x(1)^2+x(3)^2);
    xdot = [x(2); -mu*x(1)/r^3 + wk(1); x(4); -mu*x(3)/r^3 + wk(2)];
end


%% LKF

function [xhat_plus, P_plus, Pyy_kp1, innov] = lkf(tvec, xstar, F, dxhat_plus_0, R, dy, H, P_plus_0, Q, omega)
    I = eye(4);

    xhat_plus(:,1) = xstar(:,1);
    dxhat_plus(:,1) = dxhat_plus_0;
    P_plus(:,:,1) = P_plus_0;
    
    for i = 2:length(tvec)
        % Time Update
        dxhat_minus(:,i) = F(:,:,i-1)*dxhat_plus(:,i-1);
        P_minus(:,:,i) = F(:,:,i-1)*P_plus(:,:,i-1)*F(:,:,i-1)' + omega*Q*omega';

        % Measurement Update
        Pyy_kp1(:,:,i) = H(:,:,i)*P_minus(:,:,i)*H(:,:,i)' + R;
        K(:,:,i) = P_minus(:,:,i)*H(:,:,i)'*inv(Pyy_kp1(:,:,i));
        innov(:,i) = dy(:,i) - H(:,:,i)*dxhat_minus(:,i);
        dxhat_plus(:,i) = dxhat_minus(:,i) + K(:,:,i)*innov(:,i);
        P_plus(:,:,i) = (I - K(:,:,i)*H(:,:,i))*P_minus(:,:,i);
        
        % Adding up deterministic and perturbations
        xhat_plus(:,i) = xstar(:,i) + dxhat_plus(:,i);
    end
end
        

