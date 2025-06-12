clear; clc; close all;

delta_t = 10;
k = 398600;
r0 = 6678;
omega_0 = sqrt(k/r0^3);
orbital_period = 2*pi * sqrt(r0^3/k);

A_bar = [0 1 0 0; omega_0^2+2*k/r0^3 0 0 2*r0*omega_0; 0 0 0 1; 0 -2*omega_0/r0 0 0];
B_bar = [0 0; 1 0; 0 0; 0 1/r0];

% A is always square matrix
n = length(A_bar);
m = size(B_bar,2);

A_bar_hat = [A_bar B_bar; zeros([n - m, n+m])];

syms t

stm = expm(A_bar_hat*t);
stmf(t) = stm;
f(t) = stm(1:4, 1:4);
g(t) = stm(1:4, 5:6);

% Problem 1, part c
out = double(stmf(10));

% Problem 4, part a
r_dot_0 = 0;
theta_0 = omega_0 * t;
theta_dot_0 = omega_0;
x_nom(t) = [r0, r_dot_0, theta_0, theta_dot_0]';
x_nom_0 = double(x_nom(0));

delta_r0 = 10;
delta_r_dot = -0.5;
delta_theta = 0;
delta_theta_dot = 2.5e-5;
delta_x_0 = [delta_r0, delta_r_dot, delta_theta, delta_theta_dot]';

u1 = 0;
u2 = 0;
delta_u = [u1, u2]';

delta_x = [delta_x_0];
x_nominal = [x_nom_0];

time = 0:10:orbital_period;

for i = 2:length(time)
    delta_x(:,i) = double(f(time(i))) * delta_x_0 + double(g(time(i))) * delta_u;
    x_nominal(:,i) = double(x_nom(time(i))) + delta_x(:,i);
end

figure()
subplot(4,1,1)
plot(time, delta_x(1,:))
legend("\delta r", 'FontSize', 13)
ylabel("Km")

subplot(4,1,2)
plot(time, delta_x(2,:))
legend("$\delta \dot{r}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Km/s")

subplot(4,1,3)
plot(time, delta_x(3,:))
legend("\delta \theta", 'FontSize', 13)
ylabel("rad")

subplot(4,1,4)
plot(time, delta_x(4,:))
legend("$\delta \dot{\theta}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("rad/s")
xlabel("Time [seconds]")

sgtitle("Perturbation States vs Time")


figure()
subplot(4,1,1)
plot(time, x_nominal(1,:))
legend("r", 'FontSize', 13)
ylabel("Km")

subplot(4,1,2)
plot(time, x_nominal(2,:))
legend("$\dot{r}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Km/s")

subplot(4,1,3)
plot(time, x_nominal(3,:))
legend("\theta", 'FontSize', 13)
ylabel("rad")

subplot(4,1,4)
plot(time, x_nominal(4,:))
legend("$\dot{\theta}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("rad/s")
xlabel("Time [seconds]")

sgtitle("Expected System Total States vs Time")

%%
clear; clc; close all;

k = 398600;
r0 = 6678;
omega_0 = sqrt(k/r0^3);
r_dot_0 = 0;
theta_0 = 0;
theta_dot_0 = omega_0;
u1 = 0;
u2 = 0;

x0 = [r0, r_dot_0, theta_0, theta_dot_0]';

orbital_period = 2*pi * sqrt(r0^3/k);
time = 0:10:orbital_period;

x_nom = @(t,x) [x(2); x(1)*x(4)^2 - k/x(1)^2 + u1; x(4); -2*x(4)*x(2)/x(1) + 1/x(1)*u2];

% Perturbations
A_bar = [0 1 0 0; omega_0^2+2*k/r0^3 0 0 2*r0*omega_0; 0 0 0 1; 0 -2*omega_0/r0 0 0];
B_bar = [0 0; 1 0; 0 0; 0 1/r0];

dx0 =  [10; -0.5; 0; 2.5e-5];
x_perturb = @(t, dx) [dx(2); (omega_0^2+2*k/r0^3)*dx(1) + 2*r0*omega_0*dx(4); dx(4); -2*omega_0/r0*dx(2)];

[t_nom, x_nominal] = ode45(x_nom, time, x0);
[t_perturb, x_perturb] = ode45(x_perturb, time, dx0);

x = x_nominal + x_perturb;


figure()
subplot(4,1,1)
plot(time, x_nominal(:,1))
legend("r", 'FontSize', 13)
ylabel("Km")

subplot(4,1,2)
plot(time, x_nominal(:,2))
legend("$\dot{r}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Km/s")

subplot(4,1,3)
plot(time, x_nominal(:,3))
legend("\theta", 'FontSize', 13)
ylabel("rad")

subplot(4,1,4)
plot(time, x_nominal(:,4))
legend("$\dot{\theta}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("rad/s")
xlabel("Time [seconds]")

sgtitle("Nominal System States vs Time")


figure()
subplot(4,1,1)
plot(time, x_perturb(:,1))
legend("r", 'FontSize', 13)
ylabel("Km")

subplot(4,1,2)
plot(time, x_perturb(:,2))
legend("$\dot{r}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Km/s")

subplot(4,1,3)
plot(time, x_perturb(:,3))
legend("\theta", 'FontSize', 13)
ylabel("rad")

subplot(4,1,4)
plot(time, x_perturb(:,4))
legend("$\dot{\theta}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("rad/s")
xlabel("Time [seconds]")

sgtitle("Actual Perturbation States vs Time")


figure()
subplot(4,1,1)
plot(time, x(:,1))
legend("r", 'FontSize', 13)
ylabel("Km")

subplot(4,1,2)
plot(time, x(:,2))
legend("$\dot{r}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Km/s")

subplot(4,1,3)
plot(time, x(:,3))
legend("\theta", 'FontSize', 13)
ylabel("rad")

subplot(4,1,4)
plot(time, x(:,4))
legend("$\dot{\theta}$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("rad/s")
xlabel("Time [seconds]")

sgtitle("Actual Total System States vs Time")
