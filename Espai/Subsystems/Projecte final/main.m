%% Init
close all; clearvars

%% User parameters
animation = false;
do_plots = true;
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

%% Section 1
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

if animation
    generate_Animation_GPT(r_POI_eci,r_sat_eci,R_E,t,n_steps,w_E,w_sat)
end

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

%     figure
%     title('Representation')
%     plot_earth(R_E)
%     plot3(r_POI_eci(1,1),r_POI_eci(2,1),r_POI_eci(3,1), 'ro', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
%     plot3(r_sat_eci(1,1),r_sat_eci(2,1),r_sat_eci(3,1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
end

%% === 2. Attitude dynamics and control ===
% Este bloque implementa 2.1..2.4 usando cuaterniones.
% Asume variables definidas en main.m:
%   r_sat_eci   (3xN), r_POI_eci (3xN), t (1xN), w_sat (3x1 o 3xN opcional)
% Guarda resultados en arrays para posteriores gráficas.

% ---- Parameters (ajusta si lo necesitas) ----
I = diag([100, 120, 80]);    % Inercia (kg*m^2) - ajusta según tu satélite
invI = inv(I);
% Gains PDI (puedes convertir en vectores por eje)
Kp = 2.0;
Kd = 1.0;
Ki = 0.01;

% Sensor axis (en BRF) que queremos alinear con sat->POI (por ejemplo eje +Z)
b_sensor = [0; 0; 1];

% Inicial conditions (puedes cambiar q0 y omega0)
q = [1; 0; 0; 0];        % quaternion (scalar-first) que transforma body->inertial
% Si tienes w_sat inicial en tu main.m, úsalo; en caso contrario:
if exist('w_sat','var') && numel(w_sat)==3
    omega = w_sat(:);  % body angular velocity inicial (rad/s)
elseif exist('w_sat','var') && size(w_sat,2)>=1
    omega = w_sat(:,1);
else
    omega = [0;0;0.05];
end

% Preparar arrays de salida
N = length(t);
q_hist = zeros(4,N);
omega_hist = zeros(3,N);
tau_hist = zeros(3,N);
err_hist = zeros(3,N);

int_err = zeros(3,1); % integral term init

% dt (suponemos muestreo uniforme)
if length(t)>=2
    dt = t(2)-t(1);
else
    dt = 1;
end

% Run integration loop (simple Euler explícito). Si prefieres, usa RK4/ode45.
for k = 1:N
    % posiciones
    rs = r_sat_eci(:,k);
    rp = r_POI_eci(:,k);

    % ---- 2.1 Cuaternión deseado (orientar b_sensor hacia v_des) ----
    v_des_eci = rp - rs;
    v_des_eci = v_des_eci / norm(v_des_eci);

    % Convertir b_sensor (body) a ECI con q: v_sensor_eci = R(q) * b_sensor
    q0 = q(1); qv = q(2:4);
    R_q = (2*q0^2 - 1)*eye(3) + 2*(qv*qv.') + 2*q0*skew(qv);
    v_sensor_eci = R_q * b_sensor;

    % quaternion que rota v_sensor_eci -> v_des_eci (scalar-first)
    qd = quatFromTwoVectors(v_sensor_eci, v_des_eci);
    qd = quatNormalize(qd);

    % ---- 2.2 Error de actitud ----
    qinv = quatConj(q);          % q^{-1}
    q_err = quatMultiply(qd, qinv);
    q_err = quatNormalize(q_err);
    % for control, take vector part; ensure scalar positive to use shortest rotation
    if q_err(1) < 0
        q_err = -q_err;
    end
    err_vec = q_err(2:4); % 3x1

    % ---- 2.3 Control PDI -> torque ----
    int_err = int_err + err_vec * dt;
    der_err = -omega; % heuristic: reduce angular velocity to zero
    tau_ctrl = Kp * err_vec + Ki * int_err + Kd * der_err;

    % ---- 2.3 Dinámica rígida ----
    omega_dot = invI * ( tau_ctrl - cross(omega, I*omega) );

    % ---- 2.4 Cinemática de cuaterniones ----
    qdot = 0.5 * OmegaMatrix(omega) * q;

    % Integración explícita (Euler)
    omega = omega + omega_dot * dt;
    q = q + qdot * dt;
    q = quatNormalize(q);

    % Guardar historiales
    q_hist(:,k) = q;
    omega_hist(:,k) = omega;
    tau_hist(:,k) = tau_ctrl;
    err_hist(:,k) = err_vec;
end

% Almacenamos resultados para plot posterior (si do_plots == true)
if exist('do_plots','var') && do_plots
    figure;
    subplot(3,1,1); plot(t, q_hist(1,:), t, q_hist(2,:), t, q_hist(3,:), t, q_hist(4,:));
    legend('q0','q1','q2','q3'); title('Quaternion history');

    subplot(3,1,2); plot(t, omega_hist(1,:), t, omega_hist(2,:), t, omega_hist(3,:));
    legend('\omega_x','\omega_y','\omega_z'); title('Angular velocity');

    subplot(3,1,3); plot(t, vecnorm(err_hist));
    title('Attitude error norm'); xlabel('t (s)');
end

% Final: exportar variables al espacio de trabajo (opcional)
results.q_hist = q_hist;
results.omega_hist = omega_hist;
results.tau_hist = tau_hist;
results.err_hist = err_hist;
assignin('base','ADCS_results', results);

