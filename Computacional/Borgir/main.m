% BURGERS_GALERKIN_SOLVER.m
% Solves the 1D Viscous Burgers' Equation in Fourier space using the Galerkin method.
% Includes Direct Numerical Simulation (DNS) and Large Eddy Simulation (LES).

clear all; close all; clc;

% =========================================================================
% 1. PARAMETROS DE CONFIGURACIÓN
% =========================================================================
N_list = [20, 40];      % Números de modos de truncamiento N (ejemplos: N=20 y N=40)
Re = 40;                % Número de Reynolds
nu = 1/Re;              % Viscosidad molecular
tmax = 1;               % Tiempo total de simulación
CFL = 0.025;            % Coeficiente de estabilidad (CFL)

% Parámetros del modelo LES (Viscosidad de Remolino Espectral)
CK = 0.05;              % Coeficiente de viscosidad turbulenta
m = 2;                  % Exponente de la dependencia del número de onda |k|

% =========================================================================
% 2. INICIALIZACIÓN Y BUCLE DE SIMULACIÓN
% =========================================================================
figure(1); hold on;
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [1 10 100]); set(gca, 'YTick', [1e-6 1e-4 1e-2 1 1e2]);
grid on;

energy_total_dns = {};  % Almacena energía total DNS para Figura 2
time_all = {};          % Almacena vectores de tiempo

for idx = 1:length(N_list)
    N = N_list(idx);

    % Cálculo de paso de tiempo (Euler) y número de pasos
    dt = CFL / (nu * N^2);
    nt = round(tmax/dt);
    
    % Vectores de modos k: [-N, ..., N]
    k_vec = -N:N;
    K = length(k_vec);
    pk0 = N + 1; % Índice del modo k=0
    
    % Inicialización de la solución espectral |u_k| ~ 1/|k|
    uk_dns = zeros(1, K);
    for j = 1:N
        uk_dns(pk0 + j) = 1/j;
        uk_dns(pk0 - j) = conj(uk_dns(pk0 + j));
    end
    uk_dns(pk0) = 0; % Componente media nula
    
    uk_les = uk_dns; % Solución LES inicializada

    % Forzamiento: Fija el modo k=1 (condición de contorno/gran escala)
    uk1_const = uk_dns(pk0 + 1); 

    % Viscosidad de remolino espectral LES: nu_t(k) = CK * (|k|/N)^m
    nu_t = CK * (abs(k_vec) / N).^m;
    
    energy_total_dns{idx} = zeros(1, nt); 
    
    % --- Bucle de Tiempo (Euler Explícito) ---
    for t = 1:nt
        
        % 1. CÁLCULO DEL TÉRMINO NO LINEAL C(k) (Método Galerkin: Convolución Triádica)
        C_dns = compute_convolution(uk_dns, k_vec, K);
        C_les = compute_convolution(uk_les, k_vec, K);

        % 2. CÁLCULO DEL TÉRMINO DISIPATIVO D(k)
        
        % Disipación DNS: D_dns = -nu * k^2 * u_k
        D_dns = -(k_vec.^2) * nu .* uk_dns;
        
        % Disipación LES: D_les = -(nu + nu_t) * k^2 * u_k
        nu_eff_les = nu + nu_t;
        D_les = -(k_vec.^2) .* nu_eff_les .* uk_les;
        
        % 3. INTEGRACIÓN TEMPORAL: du_k/dt = -C(k) + D(k)
        
        % Integración DNS
        uk_dns = uk_dns + dt * (-C_dns + D_dns);

        % Integración LES
        uk_les = uk_les + dt * (-C_les + D_les);

        % Mantenimiento del modo k=1 (forzamiento)
        uk_dns(pk0 + 1) = uk1_const;
        uk_dns(pk0 - 1) = conj(uk1_const);
        uk_les(pk0 + 1) = uk1_const;
        uk_les(pk0 - 1) = conj(uk1_const);
        
        % Almacenamiento de la energía total (Solo DNS)
        energy_total_dns{idx}(t) = sum(abs(uk_dns).^2);

        % Mostrar progreso
        if mod(t, round(nt/10)) == 0
            fprintf('N=%d --> Progreso %.0f%%\n', N, 100*t/nt);
        end
    end

    % 4. RESULTADOS AL TIEMPO FINAL: ESPECTRO DE ENERGÍA E(k)
    
    % Extraer solo modos positivos (sin k=0)
    k_pos = k_vec(pk0+1:end);
    E_k_dns = abs(uk_dns(pk0+1:end)).^2;
    E_k_les = abs(uk_les(pk0+1:end)).^2;

    % Graficar espectros en Figura 1
    loglog(k_pos, E_k_dns, 'o-', 'LineWidth', 1.5, 'DisplayName', ['DNS N=' num2str(N)]);
    loglog(k_pos, E_k_les, 's-', 'LineWidth', 1.5, 'DisplayName', ['LES N=' num2str(N) ' (LES)']);

    % Almacenar tiempos
    time_all{idx} = dt * (0:nt-1);
end

% =========================================================================
% 3. GRÁFICA DEL ESPECTRO DE ENERGÍA (FIGURA 1)
% =========================================================================

k_ref = k_pos(1);
ref = (k_pos / k_ref).^(-2) * E_k_dns(1); % Línea de referencia k^-2 ajustada
loglog(k_pos, ref, '--k', 'LineWidth', 1.2, 'DisplayName', 'Pendiente -2 (Referencia)');

xlabel('Modo k');
ylabel('E_k');
title('Comparación de Espectros de Energía (DNS vs LES), Re=40');
legend('Location', 'southwest');
hold off;

% =========================================================================
% 4. GRÁFICA DE EVOLUCIÓN DE LA ENERGÍA TOTAL (FIGURA 2)
% =========================================================================
figure(2);
hold on;
for idx = 1:length(N_list)
    plot(time_all{idx}, energy_total_dns{idx}, 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_list(idx))]);
end
xlabel('Tiempo (t)');
ylabel('Energía Total');
title('Evolución de la Energía Total (DNS) para diferentes N');
grid on;
legend('Location', 'best');
hold off;


% =========================================================================
% 5. FUNCIÓN AUXILIAR DE CONVOLUCIÓN GALERKIN
% =========================================================================

function C_k = compute_convolution(uk, k_vec, K)
    % Calcula el término no lineal C(k) = sum_{p+q=k} u_p * (i*q) * u_q.
    % Representa el término convectivo F(u*du/dx) en el método Galerkin.
    
    C_k = zeros(1, K);
    for j = 1:K % Modo de salida k
        k = k_vec(j);
        sum_triads = 0;
        
        for p_idx = 1:K % Modo de entrada p
            p = k_vec(p_idx);
            q = k - p; % Modo q requerido para cerrar la tríada p + q = k
            q_idx = find(k_vec == q); % Busca el índice de q
            
            % Solo se suman tríadas válidas si q está en el dominio [-N, N]
            if ~isempty(q_idx)
                up = uk(p_idx);
                uq = uk(q_idx);
                % La derivada d/dx(u_q) introduce el factor i*q
                sum_triads = sum_triads + up * (1i * q) * uq; 
            end
        end
        C_k(j) = sum_triads;
    end
end