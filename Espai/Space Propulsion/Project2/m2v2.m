% M6_HET_Solver_Fixed_Stable_Plotting.m - Script que resuelve el sistema 1D de Hall Effect Thruster (HET)

% El código ha sido corregido para evitar inestabilidades numéricas (versión anterior)
% y para solucionar los warnings en los plots (versión actual).

disp('Iniciando el cálculo del HET (Hall Effect Thruster)...');

% ----------------------------------------------------
% 1. CONSTANTES FÍSICAS Y NO-DIMENSIONALES
% ----------------------------------------------------

% Constantes Físicas 
e = 1.6022e-19; % C (Carga del electrón)
m_e = 9.1034e-31; % kg (Masa del electrón)
m_i = 2.1802e-25; % kg (Masa del ion Xenon)
L = 0.03; % m (Longitud del canal)
V_T = 300; % V (Voltaje total)
V_i = 12.1; % V (Potencial de ionización de Xenon)

% Constantes de Referencia (T* = e*Vi)
T_star = e * V_i; % J
v_star = sqrt(T_star / m_i); % m/s

% Constantes No-dimensionales
M_hat2 = m_e / m_i; % M^2 (Relación de masas)
E_i_hat = 1; % E_i/T*
E_i_prime_hat = 2.5; % E_i'/T*
Q_w_hat = 0.2; % Pérdidas de pared
sigma_0 = 3.6e-20; % m^2
sigma_en = 27e-20; % m^2
sigma_en_hat = sigma_en / sigma_0; % sigma_en/sigma_0 = 7.5
sigma_ia_hat = 1; % Asumiendo sigma_ia = sigma_0 (Intercambio de carga/Accele.)

% Parámetros de la geometría y flujo (No utilizados en las ODEs actuales)
D_ch = 0.075; % m (Diámetro del canal)
W = 0.015; % m (Ancho del canal)
m_dot_T = 5e-6; % kg/s (Flujo másico total de Xenon)
eta_m = 0.95; % Eficiencia de utilización de propelente

% Perfil de Campo Magnético B(z_hat)
B_max = 0.02; % T
z_max_hat = 0.5;
w_hat = 0.3;
B_hat_fun = @(z_hat) B_max * exp(-(z_hat - z_max_hat).^2 / (2 * w_hat^2)); 

% Término de Potencial (phi_prime_hat)
% Gradiente de potencial lineal no-dimensional: d(e*phi/T*)/dz_hat
phi_prime_hat = -(V_T) / V_i; % Es negativo, -V_T/V_i.

% Condiciones Iniciales (z_hat=0, ánodo)
v_i_in_hat = 0.05;
T_e_in_hat = 1.25;
y0 = [v_i_in_hat; T_e_in_hat];

% Rango de Integración (z_hat de 0 a 1)
z_span = [0 1];

% ----------------------------------------------------
% 2. RESOLUCIÓN DE LAS ODES CON ODE45
% ----------------------------------------------------
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Se usa una función anónima para pasar las constantes a hall_ode_fun
fun_ode = @(z_hat, y) hall_ode_fun(z_hat, y, M_hat2, E_i_hat, E_i_prime_hat, ...
    Q_w_hat, sigma_en_hat, sigma_ia_hat, phi_prime_hat);

[z_out, Y] = ode45(fun_ode, z_span, y0, options);

% ----------------------------------------------------
% 3. EXTRACCIÓN DE RESULTADOS Y CÁLCULOS FINALES
% ----------------------------------------------------
v_i = Y(:, 1);
T_e = Y(:, 2);
v_i_out_hat = v_i(end);

% Empuje (F) [F = eta_m * m_dot_T * v_i_out]
v_i_out = v_i_out_hat * v_star;
F = eta_m * m_dot_T * v_i_out;

% Eficiencia Anódica (eta_A) [eta_A = eta_m * eta_E]
% eta_E = (0.5 * m_i * v_i_out^2) / (e * V_T)
eta_E = (0.5 * m_i * v_i_out^2) / (e * V_T);
eta_A = eta_m * eta_E;

% ----------------------------------------------------
% 4. RESULTADOS DE SALIDA
% ----------------------------------------------------
disp('----------------------------------------------------');
disp('RESULTADOS FINALES:');
disp(['Empuje (Thrust, F): ', num2str(F, '%.4e'), ' N']);
disp(['Eficiencia Anódica (Anode Efficiency, eta_A): ', num2str(eta_A * 100, '%.2f'), ' %']);
disp(['Velocidad de Salida Iónica (v_i_out): ', num2str(v_i_out, '%.2f'), ' m/s']);

% Para visualizar la evolución espacial:
figure;
subplot(2,1,1);
plot(z_out, v_i, 'LineWidth', 2);
% CORRECCIÓN DE PLOTTING: Usar Interpreter 'latex' y delimitadores '$'
title('Evolución de la Velocidad Iónica No-dimensional ($\hat{v}_i$)', 'Interpreter', 'latex');
xlabel('Posición No-dimensional ($\hat{z}$)', 'Interpreter', 'latex');
ylabel('$\hat{v}_i$', 'Interpreter', 'latex');
grid on;

subplot(2,1,2);
plot(z_out, T_e, 'LineWidth', 2);
% CORRECCIÓN DE PLOTTING: Usar Interpreter 'latex' y delimitadores '$'
title('Evolución de la Temperatura Electrónica No-dimensional ($\hat{T}_e$)', 'Interpreter', 'latex');
xlabel('Posición No-dimensional ($\hat{z}$)', 'Interpreter', 'latex');
ylabel('$\hat{T}_e$', 'Interpreter', 'latex');
grid on;


% ----------------------------------------------------
% DEFINICIÓN DE LA FUNCIÓN DE ODES
% ----------------------------------------------------

function dy_dz = hall_ode_fun(z_hat, y, M_hat2, E_i_hat, E_i_prime_hat, Q_w_hat, sigma_en_hat, sigma_ia_hat, phi_prime_hat)
    % hall_ode_fun - Define el sistema de ODEs del HET.

    v_i_hat = y(1);
    T_e_hat = y(2);
    
    % ----------------------------------------------------
    % ASUNCIONES Y TÉRMINOS AUXILIARES
    % ----------------------------------------------------
    
    % ASUNCIÓN 1: Densidad de neutrales (n_n_hat) constante y unitaria 
    n_n_hat = 1; 
    
    % ASUNCIÓN 2 (CORREGIDA): Sección transversal de ionización (sigma_iz_hat).
    % Modelo Lotz-like suavizado para evitar inestabilidades numéricas.
    C_Lotz = 4.0; 
    E_i_hat_eff = 1.0; % Umbral de ionización (E_i_hat)

    if T_e_hat > E_i_hat_eff 
        % Función: C * log(Te/Ei) / (Te * Ei) - Captura el umbral y el pico.
        sigma_iz_hat = C_Lotz * (log(T_e_hat / E_i_hat_eff) / (T_e_hat * E_i_hat_eff));
        sigma_iz_hat = max(0, sigma_iz_hat); 
    else
        sigma_iz_hat = 0;
    end
    
    % Velocidad térmica del electrón no-dimensional
    v_e_th_hat = sqrt(8 * T_e_hat / (pi * M_hat2));
    
    % Frecuencias de Colisión No-dimensionales 'v_prime'
    v_iz_prime_hat = n_n_hat * sigma_iz_hat * v_e_th_hat; % Ionización
    v_ia_prime_hat = n_n_hat * sigma_ia_hat * v_e_th_hat; % Aceleración/Intercambio de carga
    v_en_prime_hat = n_n_hat * sigma_en_hat * v_e_th_hat; % Momentum e-neutral
    
    % Pérdidas de Energía del Electrón No-dimensionales (E_ix_prime_hat)
    E_iz_prime_hat = v_iz_prime_hat * E_i_hat; % Pérdida por Ionización
    E_ia_prime_hat = v_ia_prime_hat * E_i_prime_hat; % Pérdida por Excitación (E_i')
    
    % E_el_prime_hat (Pérdida Elástica)
    E_el_prime_hat = (2 * M_hat2 / 3) * v_en_prime_hat * T_e_hat; 
    
    E_ix_prime_hat = E_iz_prime_hat + E_ia_prime_hat + E_el_prime_hat + Q_w_hat; % E_rad' = 0
    
    % ----------------------------------------------------
    % DEFINICIÓN DE LAS ODES (Sistema: dY/dz = [f1; f2])
    % ----------------------------------------------------
    
    % 1. d(v_i_hat)/d(z_hat)
    f1_num = (-phi_prime_hat) - (v_i_hat^2 * (v_iz_prime_hat + v_ia_prime_hat));
    
    if v_i_hat <= 1e-12 
        dv_i_dz = 0; % Estabilidad numérica: evitar división por cero cerca del ánodo
    else
        dv_i_dz = (1 / v_i_hat) * f1_num;
    end
    
    % 2. d(T_e_hat)/d(z_hat)
    f2_num_1 = (1 / T_e_hat) * (v_i_hat^2 * v_en_prime_hat / M_hat2) * T_e_hat; % Termino de aceleración J_e E
    f2_num = f2_num_1 - E_ix_prime_hat;
    
    f2_den_orig = v_iz_prime_hat + (1 / T_e_hat) * (v_ia_prime_hat / M_hat2 + v_en_prime_hat);

    if f2_den_orig <= 1e-12 || T_e_hat <= 1e-12 
        dT_e_dz = 0; 
    else
        dT_e_dz = (1 / T_e_hat) * f2_num / f2_den_orig;
    end
    
    dy_dz = [dv_i_dz; dT_e_dz];
end