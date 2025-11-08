%% Init
close all; clearvars

%% User parameters
animation = true;
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

% Ensamblaje y cálculos de elementos del satélite
Sat = [];
Sat = [Sat, Body('cube', 400, [0,0,0], 2, eye(3))];
Sat = [Sat, Body('box', 50, [3,0,0],[4,2,0.1],eye(3))];
Sat = [Sat, Body('box', 50, [-3,0,0],[4,2,0.1],eye(3))];

[r_cdg_sat_local, I_sat] = assemblySatellite(Sat);

% Parámetros orbitales
h = 700e3;
r_orbit = R_E + h;

T_orbit = 2*pi*sqrt(r_orbit^3/mu_E);
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

%% Section 2: Attitude dynamics and control

% Calculamos solo una vez la inversa
invI_sat = inv(I_sat);

% Ganancias iniciales de control (ajustar según respuesta)
Kp = 0.003;
Ki = 0;
Kd = 5;

% Eje del sensor (en el sistema de cuerpo)
b_sensor = [0; 0; 1];

% Condiciones iniciales
q = [1; 0; 0; 0];           % quaternion body->inertial
omega = w_sat(:);            % usa la velocidad inicial definida en SI

% Prealocación
q_hist = zeros(4,n_steps);
omega_hist = zeros(3,n_steps);
tau_hist = zeros(3,n_steps);
err_ang_hist = zeros(3,n_steps);
int_err = zeros(3,1);

for k = 1:n_steps
    rs = r_sat_eci(:,k);
    rp = r_POI_eci(:,k);

    % 2.1 Cuaternión deseado
    v_des_eci = rp - rs;
    v_des_eci = v_des_eci / norm(v_des_eci);

    % Vector del sensor en ECI
    q0 = q(1);
    qv = q(2:4);
    R_q = (2*q0^2 - 1)*eye(3) + 2*(qv*qv.') + 2*q0*skew(qv);
    v_sensor_eci = R_q * b_sensor;

    qd = quatFromTwoVectors(v_sensor_eci, v_des_eci);
    qd = quatNormalize(qd);

    % 2.2 Cuaternión de error
    qinv = quatConj(q);
    q_err = quatMultiply(qd, qinv);
    q_err = quatNormalize(q_err);
    if q_err(1) < 0
        q_err = -q_err; % asegurar rotación mínima
    end

    % Convertimos a axis-angle (axis * angle)
    q0_err = max(min(q_err(1),1),-1);
    vec_err = q_err(2:4);
    angle = 2 * acos(q0_err);
    if abs(angle) < 1e-8 % Si es muy pequeño evitamos hacer el tonto y que se disparen los numeros
        err_ang_vec = [0;0;0];
    else
        s = sin(angle/2);
        if abs(s) < 1e-8
            axis = vec_err / norm(vec_err);
        else
            axis = vec_err / s;
        end
        axis = axis / norm(axis);
        err_ang_vec = axis * angle;
    end
    err_ang_hist(:,k) = err_ang_vec;

    % 2.3 Control PDI
    int_err = int_err + err_ang_vec * dt;
    der_err = -omega;
    tau_ctrl = Kp * err_ang_vec + Ki * int_err + Kd * der_err;

    % 2.4 Dinámica y cinemática
    omega_dot = invI_sat * ( tau_ctrl - cross(omega, I_sat*omega) );
    qdot = 0.5 * OmegaMatrix(omega) * q;

    % Integración explícita (Euler)
    omega = omega + omega_dot * dt;
    q = q + qdot * dt;
    q = quatNormalize(q);

    % Guardar historiales
    q_hist(:,k) = q;
    omega_hist(:,k) = omega;
    tau_hist(:,k) = tau_ctrl;
end

% He probado de hacer un plot con el movimiento final del satélite pero no
% es demasiado bonito, creo que la función no está preparada para el
% formato de omega_hist
%generate_Animation_GPT(r_POI_eci,r_sat_eci,R_E,t,n_steps,w_E,omega_hist)

% Resultados
results.q_hist = q_hist;
results.omega_hist = omega_hist;
results.tau_hist = tau_hist;
results.err_ang_hist = err_ang_hist;
assignin('base','ADCS_results', results);

% --- Plots finales ---
if do_plots
    % Vector error (norma)
    figure('Name','Error vector (axis * angle)');
    plot(t, vecnorm(err_ang_hist));
    grid on; xlabel('t (s)'); ylabel('||axis*angle|| (rad)');
    title('Norma del vector de error (axis*angle)');

    % Velocidad angular
    figure('Name','Angular velocity');
    plot(t, omega_hist(1,:), t, omega_hist(2,:), t, omega_hist(3,:),'LineWidth',1.2);
    legend('\omega_x','\omega_y','\omega_z'); grid on;
    xlabel('t (s)'); ylabel('\omega (rad/s)');
    title('Velocidad angular del satélite');

    % Par de control
    figure('Name','Control torques');
    plot(t, tau_hist(1,:), t, tau_hist(2,:), t, tau_hist(3,:),'LineWidth',1.2);
    legend('\tau_x','\tau_y','\tau_z'); grid on;
    xlabel('t (s)'); ylabel('\tau (N·m)');
    title('Componentes del par de control');
end

