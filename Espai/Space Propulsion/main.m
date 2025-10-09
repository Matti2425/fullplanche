close all
clearvars
clc
clear

mu = 3.986e5;
g0 = 9.81;  %m/s^2
Isp = 2000; %s
F = 10*0.001;      %N
alpha = 0;  %deg
alpha = deg2rad(alpha);

options = odeset('AbsTol', 1e-9, 'RelTol',1e-9);

% Init initial params
ri = (6371+200);
Vtheta_i = 0.8* sqrt(mu/ri);
mi = 800;

xi = [ri, 0, 0, Vtheta_i, mi];

tspan = linspace(0,1e6,1e5+1);
[t,x] = ode45(@(t,x) ode(t,x,mu,g0,Isp,F,alpha),tspan,xi, options);

xc= x(:,1) .* cos(x(:,2));
yc= x(:,1) .* sin(x(:,2));

plot(xc,yc)
axis equal