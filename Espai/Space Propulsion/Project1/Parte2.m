% MissionToMars_Part2_Corrected.m
% Código corregido con cálculos de tiempo de vuelo

clear; close all; clc;

%% 1. DATOS Y CONSTANTES
mu_sun = 1.32712440018e20;
mu_earth = 3.986004418e14;
mu_mars = 4.282837e13;
AU = 149597870e3;
R_earth = 6371e3;
R_mars = 3389.5e3;
a_earth = 1.00000011 * AU;
a_mars = 1.523679 * AU;
m0_earth = 25000; % kg
Isp = 1000; % s
g0 = 9.80665; % m/s^2
c_e = Isp * g0;

%% 2. RESULTADOS PREVIOS (de Parte 1)
rSOI_Earth = a_earth * (5.9722e24 / 1.98847e30)^(2/5);
rSOI_Mars = a_mars * (6.4171e23 / 1.98847e30)^(2/5);

% Hohmann Transfer Delta-V
a_t = (a_earth + a_mars)/2;
v_dep_helio = sqrt(mu_sun * (2/a_earth - 1/a_t));
v_arr_helio = sqrt(mu_sun * (2/a_mars - 1/a_t));
v_earth = sqrt(mu_sun / a_earth);
v_mars = sqrt(mu_sun / a_mars);
dV1_target = abs(v_dep_helio - v_earth);
dV2_target = abs(v_mars - v_arr_helio);

fprintf('--- OBJETIVOS ---\n');
fprintf('V_inf Tierra (dV1): %.2f m/s\n', dV1_target);
fprintf('V_inf Marte (dV2) : %.2f m/s\n\n', dV2_target);

%% 3. CUESTIÓN 4: ESCAPE DE LA TIERRA
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
    error_v = v_final - dV1_target;
    
    fprintf('Iter %d: F=%.1f N -> V=%.2f m/s (error=%.2f)\n', ...
            iter, F_try, v_final, error_v);
    
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

fprintf('>>> Q4: F = %.2f N, Masa Final = %.2f kg\n', F_opt_Earth, results_Q4.m_final);
fprintf('>>> Tiempo de vuelo: %.2f días\n', time_Q4_days);
fprintf('>>> Propelente usado: %.2f kg\n\n', propellant_Q4);

%% 4. CUESTIÓN 5: CAPTURA EN MARTE
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
        err_v = v_end - dV2_target;
        
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
    fprintf('Iter Masa %d: m_dry = %.1f -> m_SOI = %.1f (Target: %.1f) | F = %.1f N\n', ...
            i_mass, m_dry_guess, m_soi_reached, target_mass_SOI, F_M);
    
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

fprintf('>>> Q5: Empuje = %.2f N, Masa Orbita Marte = %.2f kg\n', F_opt_Mars, m_dry_final);
fprintf('>>> Tiempo de vuelo: %.2f días\n', time_Q5_days);
fprintf('>>> Propelente usado: %.2f kg\n\n', propellant_Q5);

%% 5. CÁLCULO DE ÓRBITA Y TIEMPOS TOTALES

% Tiempo de transferencia heliocéntrica (de Parte 1)
TOF_hohmann = pi * sqrt(a_t^3 / mu_sun);
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

%% 6. GRÁFICOS MEJORADOS

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

% Gráfico de evolución de masa
figure('Name','Evolución de Masa','Color','w');
subplot(2,1,1);
plot(results_Q4.t/86400, results_Q4.y(:,5), 'b-', 'LineWidth', 2);
xlabel('Tiempo (días)'); ylabel('Masa (kg)');
title('Fase Tierra: Evolución de Masa');
grid on;

subplot(2,1,2);
plot(results_Q5.t/86400, results_Q5.y(:,5), 'r-', 'LineWidth', 2);
xlabel('Tiempo (días)'); ylabel('Masa (kg)');
title('Fase Marte: Evolución de Masa');
grid on;

% Gráfico de velocidades
figure('Name','Evolución de Velocidad','Color','w');
subplot(2,1,1);
v_total_Q4 = sqrt(results_Q4.y(:,3).^2 + results_Q4.y(:,4).^2);
plot(results_Q4.t/86400, v_total_Q4/1e3, 'b-', 'LineWidth', 2);
hold on;
yline(dV1_target/1e3, 'r--', 'LineWidth', 1, 'Label', 'Objetivo');
xlabel('Tiempo (días)'); ylabel('Velocidad (km/s)');
title('Fase Tierra: Velocidad Total');
legend('Velocidad', 'Objetivo', 'Location', 'best');
grid on;

subplot(2,1,2);
v_total_Q5 = sqrt(results_Q5.y(:,3).^2 + results_Q5.y(:,4).^2);
plot(results_Q5.t/86400, v_total_Q5/1e3, 'r-', 'LineWidth', 2);
hold on;
yline(dV2_target/1e3, 'r--', 'LineWidth', 1, 'Label', 'Objetivo');
xlabel('Tiempo (días)'); ylabel('Velocidad (km/s)');
title('Fase Marte: Velocidad Total');
legend('Velocidad', 'Objetivo', 'Location', 'best');
grid on;

%% FUNCIONES

function dydt = dynamics_2D(~, y, F, c_e, mu, mode)
    r = y(1);
    vr = y(3);
    vt = y(4);
    m = y(5);
    
    % Aceleraciones gravitatorias
    ar_grav = -mu/r^2 + vt^2/r;
    at_grav = -vr*vt/r;
    
    % Aceleración del empuje
    v_mod = sqrt(vr^2 + vt^2);
    if v_mod > 1e-3
        ar_thrust = (F/m) * (vr/v_mod);
        at_thrust = (F/m) * (vt/v_mod);
    else
        ar_thrust = 0;
        at_thrust = 0;
    end
    
    dydt = [vr;
            vt/r;
            ar_grav + ar_thrust;
            at_grav + at_thrust;
            -mode * F/c_e];
end

function [value, isterminal, direction] = event_Check(~, y, r_soi, r_planet)
    r = y(1);
    value1 = r - r_soi;
    value2 = r - (r_planet + 100e3);
    value = [value1, value2];
    isterminal = [1, 1];
    direction = [1, -1];
end