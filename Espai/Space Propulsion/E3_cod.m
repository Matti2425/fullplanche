
%% MAIN
% Init
clear; clc; close all;

% Constants
mu = 3.986e14;        % [m^3/s^2]
g0 = 9.80665;         % [m/s^2]
Re = 6378e3;          % [m]

Isp = 1800;           % [s]
F = 0.25;             % [N] (250 mN)
m0 = 800;             % [kg]

h0 = 200e3;           % [m]
hf = 16565.96e3;      % [m] 
r0 = Re + h0;
rf = Re + hf;

% Initial state
v0 = sqrt(mu / r0);
x0 = [ r0; 0; 0; v0; m0 ];   % [r, theta, rdot, r*theta_dot, m]

% Climb angle
alpha = deg2rad(0);   % deg

% ODE options with events (stop when r >= rf)
opts = odeset('RelTol',1e-9,'AbsTol',1e-9, ...
              'Events', @(t,x) events_reach_radius(t,x,rf));

% maximum allowed integration time (large)
tmax = 3e7;  % [s]

% Integrate using ode function
[t, x, te, xe, ie] = ode45(@(t,x) ode(t, x, mu, g0, Isp, F, alpha), [0 tmax], x0, opts);

% Check whether target was reached
if isempty(te)
    warning('Target radius not reached within tmax = %.2e s', tmax);
    t_transfer = t(end);
    x_final = x(end,:)';
else
    t_transfer = te(1);
    x_final = xe(1,:)';
end

% Extract results
r_hist = x(:,1);
alt_hist = r_hist - Re;
m_final = x_final(5);

propellant_used = m0 - m_final;
mass_fraction = propellant_used / m0;

% Print results with precision
fprintf('Transfer results (using user''s ode)\n');
fprintf('Climb angle alpha = %.6f deg\n', rad2deg(alpha));
fprintf('Transfer time = %.9f s  (%.9f hours)\n', t_transfer, t_transfer/3600);
fprintf('Propellant used = %.9f kg\n', propellant_used);
fprintf('Mass fraction = %.12f\n', mass_fraction);
fprintf('Final mass = %.9f kg\n', m_final);
fprintf('Final radius (r) = %.6f m  -> altitude = %.6f km\n', x_final(1), (x_final(1)-Re)/1e3);

% Plots
figure;
subplot(2,1,1);
plot(t/3600, alt_hist/1e3);
grid on; xlabel('Time [h]'); ylabel('Altitude [km]');
title(sprintf('Altitude vs time (alpha = %.2f deg)', rad2deg(alpha)));

subplot(2,1,2);
plot(t/3600, x(:,5));
grid on; xlabel('Time [h]'); ylabel('Mass [kg]');
title('Mass vs time');

%% FUNCTIONS
function [value, isterminal, direction] = events_reach_radius(t, x, rf)
    value = x(1) - rf;    % stop when r - rf = 0
    isterminal = 1;       % stop the integration
    direction = 1;        % only trigger when crossing from below to above (increasing)
end

