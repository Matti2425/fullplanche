%% Init
close all; clearvars

%% User parameters
do_plots = false;
n_steps = 600;

%% Initial data

G = 6.6743e-11; % m^3/kg/s
m_E = 5.972e24; % kg
R_E = 6371e3;   % m
c = 299792458;  % m/s
mu_E = G*m_E;

T_E = 86400;
w_E = 2*pi/T_E;

g = 9.8066;     % m/s^2

h = 700e3;
r_orbit = R_E + h;
Sat = [];
Sat = [Sat, Body('cube', 400, [0,0,0], 2, eye(3))];
Sat = [Sat, Body('box', 50, [3,0,0],[4,2,0.1],eye(3))];
Sat = [Sat, Body('box', 50, [-3,0,0],[4,2,0.1],eye(3))];

T_orbit = 2*pi^2*sqrt(r_orbit^3/mu_E);
w_orbit = 2*pi/T_orbit;

w_sat = [0,0,0.05];

% POI
lat = deg2rad(0);
lon = deg2rad(-78.467834);

% Time discretization
t = linspace(0,T_orbit,n_steps);
dt = t(2)-t(1);

% POI position
%r_POI = [R_E*cos(lat)*cosd(lon), R_E*cosd(lat)*sind(lon), R_E*sind(lat)];
r_POI_eci = [R_E*cos(lat)*cos(lon+w_E*t(:)'); R_E*cos(lat)*sin(lon+w_E*t(:)'); R_E*sin(lat)*t(:)']; % Este truco del z sirve solo porque =0

% Sat position
theta = w_orbit * t;

r_sat_eci = zeros(3, n_steps);
for i=1:n_steps
    r_sat_eci(:,i) = r_orbit * [ cos(theta(i)); sin(theta(i)); 0 ];
end

generate_Animation(r_POI_eci,r_sat_eci,R_E,t, n_steps)

if do_plots
    figure
    hold on;
    axis equal;
    grid on;
    view(3);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Satellite Geometry');
    plot_satellite(Sat, gca)

    figure
    title('Representation')
    plot_earth(R_E)
    plot3(r_POI_eci(1,1),r_POI_eci(2,1),r_POI_eci(3,1), 'ro', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
    plot3(r_sat_eci(1,1),r_sat_eci(2,1),r_sat_eci(3,1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
end