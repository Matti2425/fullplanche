% ANTONIO LUQUE LÓPEZ
% MUEA 
% ASSIGNMENT 4: BURGERS EQUATION
%--------------------------------------------------------------------------

clear all; close all; clc;

%GENERAL PARAMETERS

%%Please, don't try high N. A lot of time to compile.

N_list = [20, 40];          
Re = 40;                   
nu = 1/Re;                 
CFL = 0.025;                 % CFL number dt
tmax = 1;                    % Simulation time

% LES PARAMETERS
CK = 0.05;                  % Turbulent viscosity coefficient
m = 2;                     

% GRAPHIC CONFIGURATIONS
figure;
hold on;
set(gca, 'XScale', 'log', 'YScale', 'log')   
set(gca, 'XTick', [1 10 100])                 
set(gca, 'YTick', [1e-6 1e-4 1e-2 1 1e2])     
grid on;

% Prepare to save total energy
energy_total_all = {};  
time_all = {};            

for idx = 1:length(N_list)
    N = N_list(idx);

    dt = CFL / (nu*N^2);
    nt = round(tmax/dt);

    energy_total = zeros(1, nt);  % Initialize total energy for this N

    k_vec = -N:N;
    K = length(k_vec);
    pk0 = N+1;

    %INITIAL CONDITIONS
    uk_dns = zeros(1, K);
    for j = 1:N
        uk_dns(pk0 + j) = 1/j;
    end
    for j = 1:N
        uk_dns(pk0 - j) = conj(uk_dns(pk0 + j));
    end
    uk_dns(pk0) = 0; % Set u_0 = 0 (no mean component)
    uk_les = uk_dns; % Initialize LES solution equal to DNS

    uk1_const_dns = uk_dns(pk0 + 1);
    uk1_const_les = uk_les(pk0 + 1);

    %TURBULENT VISCOSITY FOR LES
    nu_t = CK * (abs(k_vec)/N).^m;


    %TIME LOOP
    for t = 1:nt
        %DNS
        C_dns = zeros(1, K);
        for j = 1:K
            k = k_vec(j);
            sum_triads = 0;
            for p_idx = 1:K
                p = k_vec(p_idx);
                q = k - p;
                q_idx = find(k_vec == q);
                if ~isempty(q_idx)
                    up = uk_dns(p_idx);
                    uq = uk_dns(q_idx);
                    sum_triads = sum_triads + up * 1i * q * uq;
                end
            end
            C_dns(j) = sum_triads;
        end

        D_dns = -(k_vec.^2) / Re .* uk_dns;
        duk_dt_dns = -C_dns + D_dns;
        uk_dns = uk_dns + dt * duk_dt_dns;

        uk_dns(pk0 + 1) = uk1_const_dns;
        uk_dns(pk0 - 1) = conj(uk1_const_dns);

        %LES
        C_les = zeros(1, K);
        for j = 1:K
            k = k_vec(j);
            sum_triads = 0;
            for p_idx = 1:K
                p = k_vec(p_idx);
                q = k - p;
                q_idx = find(k_vec == q);
                if ~isempty(q_idx)
                    up = uk_les(p_idx);
                    uq = uk_les(q_idx);
                    sum_triads = sum_triads + up * 1i * q * uq;
                end
            end
            C_les(j) = sum_triads;
        end

        D_les = -(k_vec.^2) .* (nu + nu_t) .* uk_les;
        duk_dt_les = -C_les + D_les;
        uk_les = uk_les + dt * duk_dt_les;

        uk_les(pk0 + 1) = uk1_const_les;
        uk_les(pk0 - 1) = conj(uk1_const_les);

        % Save total energy at each time step
        energy_total(t) = sum(abs(uk_dns).^2);

        % Progress display
        if mod(t, round(nt/100)) == 0
            fprintf('N=%d --> Progress %.1f%%\n', N, 100*t/nt);
        end
    end

    %ENERGY SPECTRUM
    k_pos = k_vec(pk0+1:end);
    E_k_dns = abs(uk_dns(pk0+1:end)).^2;
    E_k_les = abs(uk_les(pk0+1:end)).^2;

    %PLOT SPECTRUM
    loglog(k_pos, E_k_dns, 'o-', 'LineWidth', 1.5, 'DisplayName', ['DNS N=' num2str(N)]);
    loglog(k_pos, E_k_les, 's-', 'LineWidth', 1.5, 'DisplayName', ['LES N=' num2str(N)]);

    %Total energy and time vectors
    energy_total_all{idx} = energy_total;
    time_all{idx} = dt * (0:nt-1);
end

%REFERENCE SLOPE -2 
k_ref = k_pos(1);  
ref = (k_pos / k_ref).^(-2) * E_k_dns(1);
loglog(k_pos, ref, '--k', 'LineWidth', 1.2, 'DisplayName', 'Slope -2');

xlabel('Mode k');
ylabel('E_k');
title('DNS vs LES for different N values, Re=40');
legend('Location', 'southwest');
hold off;

%PLOT TOTAL ENERGY EVOLUTION FOR EACH N
figure;
hold on;
for idx = 1:length(N_list)
    plot(time_all{idx}, energy_total_all{idx}, 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_list(idx))]);
end
xlabel('Time');
ylabel('Total Energy');
title('Total energy evolution for different N values');
grid on;
legend('Location', 'best');
hold off;

