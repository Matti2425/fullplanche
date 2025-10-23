% MissionToMars_Part1.m
% MATLAB script (and local helper functions) to solve PART 1 of the project.
% - Written to be idiomatic MATLAB, modular and function-oriented.
% - Includes a scaled graphical representation of orbits and transfer.

clear; close all; clc;

%% Constants (SI)
mu_sun = 1.32712440018e20;        % m^3/s^2 (Sun)
mass_sun = 1.98847e30;           % kg
mass_earth = 5.9722e24;          % kg
mass_mars = 6.4171e23;           % kg

AU = 149597870e3;                % m
R_earth = 6371e3;                % m
R_mars  = 3389.5e3;              % m

a_earth = 1.00000011 * AU;       % m
a_mars  = 1.523679 * AU;         % m

Earth = makeBody('Earth', mass_earth, mu_sun, a_earth, R_earth);
Mars  = makeBody('Mars' , mass_mars , mu_sun, a_mars , R_mars);
Sun   = makeBody('Sun'  , mass_sun  , mu_sun, 0       , 696340e3);

%% 1) Orbital elements
fprintf('PART 1 - Q1: Orbital elements (circular, coplanar)\n');
fprintf('Earth: a = %.3e km\n', Earth.a/1e3);
fprintf('Mars : a = %.3e km\n\n', Mars.a/1e3);

%% 2) Radii of spheres of influence
rSOI_Earth = computeSOI(Earth.mass, Sun.mass, Earth.a);
rSOI_Mars  = computeSOI(Mars.mass , Sun.mass, Mars.a );

fprintf('PART 1 - Q2: Spheres of Influence\n');
fprintf('Earth SOI = %.3f km\n', rSOI_Earth/1e3);
fprintf('Mars  SOI = %.3f km\n\n', rSOI_Mars/1e3);

%% 3) Hohmann transfer
hoh = hohmannTransfer(mu_sun, Earth.a, Mars.a);
phase_rad = requiredPhaseAngle(mu_sun, Earth.a, Mars.a);

fprintf('PART 1 - Q3: Hohmann transfer results\n');
fprintf('a_t = %.3e km\n', hoh.a_t/1e3);
fprintf('TOF = %.2f days\n', hoh.t_to_flight/86400);
fprintf('dV1 = %.3f km/s, dV2 = %.3f km/s\n', hoh.dV1/1e3, hoh.dV2/1e3);
fprintf('Phase angle = %.2f deg\n\n', rad2deg(phase_rad));

%% 4) Graphical representation (scaled)
figure('Name','Earth–Mars Hohmann Transfer','Color','w'); hold on; axis equal;

% Orbits
plotOrbit(a_earth, 'b', 'Earth orbit');
plotOrbit(a_mars,  'r', 'Mars orbit');
plotOrbit(hoh.a_t, '--k', 'Transfer orbit', Earth.a, Mars.a);

% Bodies
plot(0,0,'yo','MarkerFaceColor','y','MarkerSize',12); % Sun
plot(a_earth,0,'bo','MarkerFaceColor','b');
plot(a_mars*cos(phase_rad), a_mars*sin(phase_rad),'ro','MarkerFaceColor','r');

legend('Location','bestoutside');
xlabel('x (m)'); ylabel('y (m)'); title('Hohmann Transfer: Earth to Mars');
grid on;




