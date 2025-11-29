% MissionToMars_
% MATLAB script to solve  the project.

clear; close all; clc;

%% Constants (SI)
mu_sun = 1.32712440018e20;       % m^3/s^2 (Sun)
mu_earth = 3.986004418e14;       % m^3/s^2 (Earth)
mu_mars = 4.282837e13;           % m^3/s^2 (Mars)
mass_sun = 1.98847e30;           % kg
mass_earth = 5.9722e24;          % kg
mass_mars = 6.4171e23;           % kg

AU = 149597870e3;                % m
R_earth = 6371e3;                % m
R_mars  = 3389.5e3;              % m

a_earth = 1.00000011 * AU;       % m
a_mars  = 1.523679 * AU;         % m

m0_earth = 30000;                % kg
Isp = 900;                      % s
g0 = 9.80665;                    % m/s^2
c_e = Isp * g0;

Earth = makeBody('Earth', mass_earth, mu_sun, a_earth, R_earth);
Mars  = makeBody('Mars' , mass_mars , mu_sun, a_mars , R_mars);
Sun   = makeBody('Sun'  , mass_sun  , mu_sun, 0       , 696340e3);

% CONSTANTS FOR PART 3
muM = mu_mars;
Isp_sat = 3000;                   
c_e_sat = Isp_sat * g0;           

F_sat = 0.2;                      
m0_sat = 200;                     

alpha_max = deg2rad(60);          
delta_i_total = deg2rad(30);      

T_rot_Mars = 1.02595675 * 86400;  % Mars sideral rotation period


%% ---------------Parte 1-----------------------------------------

fprintf('\n================ PART 1 ================\n');


%% CUESTIÓN 1: Orbital elements
fprintf('PART 1 - Q1: Orbital elements (circular, coplanar)\n');
fprintf('Earth: a = %.3e km\n', Earth.a/1e3);
fprintf('Mars : a = %.3e km\n\n', Mars.a/1e3);

%% CUESTIÓN 2: Radii of spheres of influence
rSOI_Earth = computeSOI(Earth.mass, Sun.mass, Earth.a);
rSOI_Mars  = computeSOI(Mars.mass , Sun.mass, Mars.a );

fprintf('PART 1 - Q2: Spheres of Influence\n');
fprintf('Earth SOI = %.3f km\n', rSOI_Earth/1e3);
fprintf('Mars  SOI = %.3f km\n\n', rSOI_Mars/1e3);

%% CUESTIÓN 3: Hohmann transfer
hoh = hohmannTransfer(mu_sun, Earth.a, Mars.a);
phase_rad = requiredPhaseAngle(mu_sun, Earth.a, Mars.a);

fprintf('PART 1 - Q3: Hohmann transfer results\n');
fprintf('a_t = %.3e km\n', hoh.a_t/1e3);
fprintf('TOF = %.2f days\n', hoh.t_to_flight/86400);
fprintf('dV1 = %.3f km/s, dV2 = %.3f km/s\n', hoh.dV1/1e3, hoh.dV2/1e3);
fprintf('Phase angle = %.2f deg\n\n', rad2deg(phase_rad));

%% ---------------Parte 2-----------------------------------------

fprintf('\n================ PART 2 ================\n');

%% CUESTIÓN 4: ESCAPE DE LA TIERRA
fprintf('--- INICIANDO Q4: ESCAPE DE LA TIERRA ---\n');

r_park_E = 20000e3;
v_park_E = sqrt(mu_earth / r_park_E);
y0_E = [r_park_E; 0; 0; v_park_E; m0_earth];

% Parámetros de búsqueda
F_opt_Earth = NaN;
F_min = 50;
F_max = 200;
max_iter = 30;
velocity_tolerance = 0.5;
converged_Q4 = false;

for iter = 1:max_iter
    F_try = (F_min + F_max) / 2;
    
    options = odeset('Events', @(t,y) event_Check(t,y,rSOI_Earth, R_earth), ...
                     'RelTol', 1e-6, 'AbsTol', 1e-9);
    
    [t_out, y_out, te, ye, ie] = ode45(@(t,y) dynamics_2D(t,y,F_try,c_e,mu_earth,1), ...
                                        [0 1e8], y0_E, options);
    
    if ~isempty(ie)
        if any(ie == 2)
            fprintf('Iter %d: F=%.1f N -> COLISIÓN (aumentar F)\n', iter, F_try);
            F_min = F_try;
            continue;
        end
    else
        fprintf('Iter %d: F=%.1f N -> NO LLEGÓ A SOI (aumentar F)\n', iter, F_try);
        F_min = F_try;
        continue;
    end
    
    v_final = sqrt(y_out(end,3)^2 + y_out(end,4)^2);
    error_v = v_final - hoh.dV1;
    
    %fprintf('Iter %d: F=%.1f N -> V=%.2f m/s (error=%.2f)\n', iter, F_try, v_final, error_v);
    
    if abs(error_v) < velocity_tolerance
        F_opt_Earth = F_try;
        converged_Q4 = true;
        break;
    elseif error_v < 0
        F_min = F_try;
    else
        F_max = F_try;
    end
end

if ~converged_Q4
    warning('No se convergió en Q4. Usando último valor F.');
    F_opt_Earth = F_try;
end

% Calcular tiempo de vuelo y propelente usado
time_Q4_days = t_out(end) / 86400;
propellant_Q4 = m0_earth - y_out(end,5);
m_dot_Q4 = F_opt_Earth / c_e;
theoretical_time_Q4 = propellant_Q4 / m_dot_Q4 / 86400;

results_Q4.F = F_opt_Earth;
results_Q4.m_final = y_out(end,5);
results_Q4.y = y_out;
results_Q4.t = t_out;
results_Q4.time_days = time_Q4_days;
results_Q4.propellant_used = propellant_Q4;

fprintf(' PART 2 - Q4: F = %.2f N, Masa Final = %.2f kg\n', F_opt_Earth, results_Q4.m_final);
fprintf(' Tiempo de vuelo: %.2f días\n', time_Q4_days);
fprintf(' Propelente usado: %.2f kg\n\n', propellant_Q4);

%% CUESTIÓN 5: CAPTURA EN MARTE
fprintf('--- INICIANDO Q5: CAPTURA EN MARTE ---\n');

r_circ_M = R_mars + 30000e3;
v_circ_M = sqrt(mu_mars / r_circ_M);

m_dry_guess = 12500;
target_mass_SOI = results_Q4.m_final;

F_opt_Mars = NaN;
m_dry_final = NaN;
converged_Q5 = false;

mass_tolerance = 10;
velocity_tolerance_Mars = 1.0;

for i_mass = 1:15
    F_low = 50;
    F_high = 200;
    F_M = NaN;
    y0_M = [r_circ_M; 0; 0; v_circ_M; m_dry_guess];
    
    for i_force = 1:25
        F_try = (F_low + F_high)/2;
        options = odeset('Events', @(t,y) event_Check(t,y,rSOI_Mars, R_mars), ...
                         'RelTol', 1e-6, 'AbsTol', 1e-9);
        
        [tm, ym, te, ye, ie] = ode45(@(t,y) dynamics_2D(t,y,F_try,c_e,mu_mars,-1), ...
                                      [0 1e8], y0_M, options);
        
        if ~isempty(ie)
            if any(ie == 2)
                F_high = F_try;
                continue;
            end
        else
            F_low = F_try;
            continue;
        end
        
        v_end = sqrt(ym(end,3)^2 + ym(end,4)^2);
        err_v = v_end - hoh.dV2;
        
        if abs(err_v) < velocity_tolerance_Mars
            F_M = F_try;
            m_soi_reached = ym(end,5);
            break;
        elseif err_v < 0
            F_low = F_try;
        else
            F_high = F_try;
        end
    end
    
    if isnan(F_M)
        m_dry_guess = m_dry_guess - 500;
        continue;
    end
    
    err_m = m_soi_reached - target_mass_SOI;
    %fprintf('Iter Masa %d: m_dry = %.1f -> m_SOI = %.1f (Target: %.1f) | F = %.1f N\n', i_mass, m_dry_guess, m_soi_reached, target_mass_SOI, F_M);
    
    if abs(err_m) < mass_tolerance
        m_dry_final = m_dry_guess;
        F_opt_Mars = F_M;
        results_Q5.t = tm;
        results_Q5.y = ym;
        converged_Q5 = true;
        break;
    else
        m_dry_guess = m_dry_guess - 0.5 * err_m;
    end
end

if ~converged_Q5
    warning('No se convergió en Q5.');
    F_opt_Mars = F_M;
    results_Q5.t = tm;
    results_Q5.y = ym;
    m_dry_final = m_dry_guess;
end

% Calcular tiempo de vuelo y propelente usado para Q5
time_Q5_days = results_Q5.t(end) / 86400;
propellant_Q5 = target_mass_SOI - m_dry_final;
m_dot_Q5 = F_opt_Mars / c_e;
theoretical_time_Q5 = propellant_Q5 / m_dot_Q5 / 86400;

fprintf('PART 2 - Q5: Empuje = %.2f N, Masa Orbita Marte = %.2f kg\n', F_opt_Mars, m_dry_final);
fprintf(' Tiempo de vuelo: %.2f días\n', time_Q5_days);
fprintf(' Propelente usado: %.2f kg\n\n', propellant_Q5);

%% RESULTADOS CÁLCULO DE ÓRBITA Y TIEMPOS TOTALES

% Tiempo de transferencia heliocéntrica (de Parte 1)
TOF_hohmann = pi * sqrt(hoh.a_t^3 / mu_sun);
TOF_hohmann_days = TOF_hohmann / 86400;

% Tiempos totales de misión
total_mission_time = time_Q4_days + TOF_hohmann_days + time_Q5_days;
total_propellant = propellant_Q4 + propellant_Q5;

fprintf('=== RESUMEN COMPLETO DE LA MISIÓN ===\n');
fprintf('FASE TIERRA:\n');
fprintf('  Empuje: %.1f N\n', results_Q4.F);
fprintf('  Masa inicial: %.1f kg\n', m0_earth);
fprintf('  Masa al salir SOI Tierra: %.1f kg\n', results_Q4.m_final);
fprintf('  Propelente usado: %.1f kg\n', propellant_Q4);
fprintf('  Tiempo de escape: %.2f días\n', time_Q4_days);

fprintf('\nTRANSFERENCIA HELIOCÉNTRICA:\n');
fprintf('  Tiempo de vuelo: %.2f días\n', TOF_hohmann_days);

fprintf('\nFASE MARTE:\n');
fprintf('  Empuje: %.1f N\n', F_opt_Mars);
fprintf('  Masa al llegar SOI Marte: %.1f kg\n', target_mass_SOI);
fprintf('  Masa seca final: %.1f kg\n', m_dry_final);
fprintf('  Propelente usado: %.1f kg\n', propellant_Q5);
fprintf('  Tiempo de captura: %.2f días\n', time_Q5_days);

fprintf('\nTOTAL MISIÓN:\n');
fprintf('  Tiempo total: %.2f días\n', total_mission_time);
fprintf('  Propelente total: %.1f kg\n', total_propellant);
fprintf('  Masa final/orbita Marte: %.1f kg\n', m_dry_final);


%% GRAFICOS

%PARTE 1
figure('Name','Earth–Mars Hohmann Transfer','Color','w'); hold on; axis equal;

% Orbits
plotOrbit(a_earth, 'b', 'Earth orbit');
plotOrbit(a_mars,  'r', 'Mars orbit');
plotOrbit(hoh.a_t, '--g', 'Transfer orbit', Earth.a, Mars.a);

% Bodies
plot(0,0,'yo','MarkerFaceColor','y','MarkerSize',12); % Sun
plot(a_earth,0,'bo','MarkerFaceColor','b');
plot(a_mars*cos(phase_rad), a_mars*sin(phase_rad),'ro','MarkerFaceColor','r');

legend('Location','bestoutside');
xlabel('x (m)'); ylabel('y (m)'); title('Hohmann Transfer: Earth to Mars');
grid on;


%PARTE 2

% Gráfico de trayectorias
figure('Name','Trayectorias Completa','Color','w', 'Position', [100 100 1200 500]);

% Escape Tierra
subplot(1,2,1);
polarplot(results_Q4.y(:,2), results_Q4.y(:,1)/1e3, 'b-', 'LineWidth', 2);
title({'Q4: Escape Tierra', sprintf('F=%.1f N, m_{final}=%.0f kg', results_Q4.F, results_Q4.m_final)});
rlim([0 rSOI_Earth/1e3 * 1.1]);

% Captura Marte
subplot(1,2,2);
polarplot(results_Q5.y(:,2), results_Q5.y(:,1)/1e3, 'r-', 'LineWidth', 2);
title({'Q5: Captura Marte', sprintf('F=%.1f N, m_{dry}=%.0f kg', F_opt_Mars, m_dry_final)});
rlim([0 rSOI_Mars/1e3 * 1.1]);
hold on;
theta_circle = linspace(0, 2*pi, 500);
polarplot(theta_circle, r_circ_M/1000 * ones(size(theta_circle)), ...
          'r--', 'LineWidth', 1.8);

legend('Trayectoria', 'Órbita', 'Location','best');

% Gráfico de evolución de masa
figure('Name','Evolución de Masa','Color','w');
subplot(2,1,1);
plot(results_Q4.t/86400, results_Q4.y(:,5), 'b-', 'LineWidth', 2);
xlabel('Tiempo (días)'); ylabel('Masa (kg)');
title('Fase Tierra: Evolución de Masa');
grid on;

subplot(2,1,2); 
% Usar la masa final de Tierra como referencia
initial_mass_earth_end = results_Q4.y(end,5);  % Masa al final de fase Tierra
initial_mass_marte = results_Q5.y(1,5);        % Masa inicial en datos de Marte

% Calcular cuánto "aumentó" incorrectamente la masa en la integración
mass_increase = results_Q5.y(:,5) - initial_mass_marte;

% Corregir: empezar desde la masa final de Tierra y restar el aumento
mass_corrected = initial_mass_earth_end - mass_increase;

plot(results_Q5.t/86400, mass_corrected, 'r-', 'LineWidth', 2);
xlabel('Tiempo (días)');
ylabel('Masa (kg)');
title('Fase Marte: Evolución de Masa');
grid on;

%{
subplot(2,1,2);
plot(results_Q5.t/86400, results_Q5.y(:,5), 'r-', 'LineWidth', 2);
xlabel('Tiempo (días)'); ylabel('Masa (kg)');
title('Fase Marte: Evolución de Masa');
grid on;
%}

% Gráfico de velocidades
figure('Name','Evolución de Velocidad','Color','w');
subplot(2,1,1);
v_total_Q4 = sqrt(results_Q4.y(:,3).^2 + results_Q4.y(:,4).^2);
plot(results_Q4.t/86400, v_total_Q4/1e3, 'b-', 'LineWidth', 2);
hold on;
yline(hoh.dV1/1e3, 'r--', 'LineWidth', 1, 'Label', 'Objetivo');
xlabel('Tiempo (días)'); ylabel('Velocidad (km/s)');
title('Fase Tierra: Velocidad Total');
legend('Velocidad', 'Objetivo', 'Location', 'best');
grid on;

subplot(2,1,2);
v_total_Q5 = sqrt(results_Q5.y(:,3).^2 + results_Q5.y(:,4).^2);
plot(results_Q5.t/86400, v_total_Q5/1e3, 'r-', 'LineWidth', 2);
hold on;
yline(hoh.dV2/1e3, 'r--', 'LineWidth', 1, 'Label', 'Objetivo');
xlabel('Tiempo (días)'); ylabel('Velocidad (km/s)');
title('Fase Marte: Velocidad Total');
legend('Velocidad', 'Objetivo', 'Location', 'best');
grid on;


% ANÁLISIS DE VELOCIDAD - Q4 (ESCAPE DE LA TIERRA)
figure('Name','Q4: Análisis de Velocidad (Tierra)','Color','w', 'Position', [100 100 800 600]);

% Datos Q4
t_days_Q4 = results_Q4.t / 86400;
Vr_Q4 = results_Q4.y(:,3); 
Vt_Q4 = results_Q4.y(:,4); 
V_tot_Q4 = sqrt(Vr_Q4.^2 + Vt_Q4.^2); % Velocidad Total

% Subplot 1: Velocidad Tangencial
subplot(2,1,1);
plot(t_days_Q4, V_tot_Q4/1e3, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Velocidad Total |V|');
hold on;
plot(t_days_Q4, Vt_Q4/1e3, 'b-', 'LineWidth', 1, 'DisplayName', 'V_{\theta} (Tangencial)');
ylabel('Velocidad [km/s]'); 
xlabel('Tiempo [días]');
title('Q4 (Tierra): Evolución de la Velocidad Total');
legend('Location','best'); 
grid on;

% Subplot 2: Velocidad Radial 
subplot(2,1,2);
plot(t_days_Q4, V_tot_Q4/1e3, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Velocidad Total |V|');
hold on;
plot(t_days_Q4, Vr_Q4/1e3, 'r-', 'LineWidth', 1, 'DisplayName', 'V_r (Radial)');
ylabel('V_r [km/s]'); 
xlabel('Tiempo [días]');
title('Q4 (Tierra): Componente Radial (Oscilaciones debido al empuje)');
legend('Location','best'); 
grid on;


% ANÁLISIS DE VELOCIDAD - Q5 (CAPTURA EN MARTE)
figure('Name','Q5: Análisis de Velocidad (Marte)','Color','w', 'Position', [950 100 800 600]);

% Datos Q5
t_days_Q5 = results_Q5.t / 86400;
Vr_Q5 = results_Q5.y(:,3); 
Vt_Q5 = results_Q5.y(:,4); 
V_tot_Q5 = sqrt(Vr_Q5.^2 + Vt_Q5.^2); % Velocidad Total

% Subplot 1: Velocidad  Tangencial
subplot(2,1,1);
plot(t_days_Q5, V_tot_Q5/1e3, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Velocidad Total |V|');
hold on;
plot(t_days_Q5, Vt_Q5/1e3, 'b-', 'LineWidth', 1, 'DisplayName', 'V_{\theta} (Tangencial)');
ylabel('Velocidad [km/s]'); 
xlabel('Tiempo de Simulación (Desde Órbita Baja hacia SOI) [días]');
title('Q5 (Marte): Evolución de la Velocidad Total');
legend('Location','best'); 
grid on;

% Subplot 2: Velocidad Radial
subplot(2,1,2);
plot(t_days_Q5, V_tot_Q5/1e3, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Velocidad Total |V|');
hold on;

plot(t_days_Q5, Vr_Q5/1e3, 'r-', 'LineWidth', 1, 'DisplayName', 'V_r (Radial)');
ylabel('Velocidad [km/s]'); 
xlabel('Tiempo de Simulación [días]');
title('Q5 (Marte): Componente Radial (Oscilaciones)');
legend('Location','best'); 
grid on;


%% ---------------Parte 3-----------------------------------------

fprintf('\n================ PART 3 ================\n');


%% CUESTIÓN 6: SYNCHRONOUS ORBIT
T_sync = T_rot_Mars;  
r_sync = (muM * (T_sync/(2*pi))^2 )^(1/3);

fprintf('\n--- Q6: SYNCHRONOUS ORBIT ---\n');
fprintf('Synchronous radius r_sync = %.6e m   (%.3f km)\n', r_sync, r_sync/1e3);

%% INITIAL CONDITIONS
i0 = deg2rad(30);
theta0 = pi/2 - i0;
phi0   = 0;

Vr0 = 0;
% Vtheta0 = 0; % V1?
% Vphi0 = sqrt(muM / r_sync);
Vtheta0 = sqrt(muM / r_sync); % V1?
Vphi0 = 0;

y0_sat = [r_sync; theta0; phi0; Vr0; Vtheta0; Vphi0; m0_sat];

fprintf('\n--- INITIAL CONDITIONS ---\n');
fprintf('r(0)     = %.6e m   (%.3f km)\n', r_sync, r_sync/1e3);
fprintf('theta(0) = %.6f rad\n', theta0);
fprintf('phi(0)   = %.6f rad\n', phi0);
fprintf('Vr(0)    = %.6f km/s\n', Vr0/1e3);
fprintf('Vtheta(0)= %.6f km/s\n', Vtheta0/1e3);
fprintf('Vphi(0)  = %.6f km/s\n', Vphi0/1e3);
fprintf('m(0)     = %.3f kg\n\n', m0_sat);

%% ============== SPHERICAL ODE ===================
function dy = satDynamics(t,y,mu,F,ce,alpha_max)
    r=y(1); th=y(2); ph=y(3);
    Vr=y(4); Vt=y(5); Vp=y(6);
    m=y(7);

    V = sqrt(Vr^2+Vt^2+Vp^2);
    alpha = alpha_max*cos(ph);

    % F_in = F*cos(alpha);
    % F_perp = F*sin(alpha);
    % 
    % Fr = F_in*(Vr/max(V,1e-12));
    % % Ftheta = 0;
    % Fphi = F_perp;

    Ftheta = F*cos(alpha);
    Fr = F*sin(alpha);
    Fphi = 0; %?  
   
    dy=zeros(7,1);

    dy(1)=Vr;
    dy(2)=Vt/r;
    dy(3)=Vp/(r*sin(th));

    dy(4)= -mu/r^2 + (Vt^2+Vp^2)/r + Fr/m;
    dy(5)= -(Vr*Vt)/r + Vp^2/(r*tan(th)) + Ftheta/m;
    dy(6)= -(Vr*Vp)/r - (Vt*Vp)/(r*tan(th)) + Fphi/m;

    dy(7)= -F/ce;
end

%% EVENT: stop when inclination = 0°
function [value,isterminal,direction] = eventInPlane(t,y)
    i = pi/2 - y(2);
    value = i;              
    isterminal = 1;
    direction = -1;
end

options_sat = odeset('Events',@eventInPlane,'RelTol',1e-9,'AbsTol',1e-12);

fprintf('Integrating until inclination reaches zero...\n');

[t_sat, y_sat] = ode45(@(t,y) satDynamics(t,y,muM,F_sat,c_e_sat,alpha_max), ...
                       [0 5e7], y0_sat, options_sat);

%% FINAL STATE
r_end = y_sat(end,1);
theta_end = y_sat(end,2);
phi_end = y_sat(end,3);
Vr_end = y_sat(end,4);
Vtheta_end = y_sat(end,5);
Vphi_end = y_sat(end,6);
m_end = y_sat(end,7);
i_end = pi/2 - theta_end;

fprintf('\n--- FINAL STATE (COMPARE WITH EXPECTED RESULTS) ---\n');
fprintf('t_end        = %.6e s  (%.6f days)\n', t_sat(end), t_sat(end)/86400);
fprintf('r(t_end)     = %.6e m   (%.6f km)\n', r_end, r_end/1e3);
fprintf('theta(t_end) = %.6f rad\n', theta_end);
fprintf('phi(t_end)   = %.6f rad\n', phi_end);
fprintf('Vr(t_end)    = %.6f km/s\n', Vr_end/1e3);
fprintf('Vtheta(t_end)= %.6e km/s\n', Vtheta_end/1e3);
fprintf('Vphi(t_end)  = %.6f km/s\n', Vphi_end/1e3);
fprintf('m(t_end)     = %.6f kg\n', m_end);
fprintf('i(t_end)     = %.6f deg\n\n', rad2deg(i_end));

%% STORE
results_P3.t = t_sat;
results_P3.y = y_sat;

%% CUESTIÓN 7: PHASE SHIFT ---
fprintf('\n--- Q7: Phase shift 20 degrees eastward ---\n');

phi_shift = deg2rad(20);
aT = F_sat / m_end;

t_phase = phi_shift / aT;

fprintf('Phase shift time = %.6f s (%.4f hours)\n', t_phase, t_phase/3600);
