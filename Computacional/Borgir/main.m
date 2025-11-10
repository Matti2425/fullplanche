% Solves the 1D viscous Burgers equation in Fourier space using Galerkin.
% Includes DNS and LES with spectral eddy-viscosity.

clear all; close all; clc;

%% 1 Configuration Parameters
% Read simulation configuration from input function
[N_list, Re, tmax, CFL, CK, m, usr] = readInputData();
nu = 1/Re; % Molecular viscosity

%% 2 Simulation Loop
if usr.doPlots
    figure(1); hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    set(gca, 'XTick', [1 10 100]); set(gca, 'YTick', [1e-6 1e-4 1e-2 1 1e2]);
    grid on;
end

% Preallocating energy and time cells
energy_total_dns = cell(1, length(N_list));
time_all = cell(1, length(N_list));

% Main loop
for idx = 1:length(N_list)
    N = N_list(idx);

    dt = CFL / (nu * N^2); % Time step based on stability condition
    nt = round(tmax/dt);   % Number of time steps

    % Wavenumber vector (Fourier modes)
    k_vec = -N:N;
    K = length(k_vec);
    pk0 = N + 1; % Index of zero mode

    % Initialize spectral solution with |u_k| appx 1/|k|
    uk_dns = zeros(1, K);
    for j = 1:N
        uk_dns(pk0 + j) = 1/j;
        uk_dns(pk0 - j) = conj(uk_dns(pk0 + j));
    end
    uk_dns(pk0) = 0; % Remove mean flow

    % LES solution initially matches DNS
    uk_les = uk_dns;

    % Store large-scale forcing (mode k=1)
    uk1_const = uk_dns(pk0 + 1);

    % Spectral eddy viscosity model for LES
    nu_t = CK * (abs(k_vec) / N).^m;

    % Preallocate energy storage
    energy_total_dns{idx} = zeros(1, nt);

    % Time integration loop
    for t = 1:nt
        % Nonlinear term (Galerkin convolution)
        C_dns = compute_convolution(uk_dns, k_vec, K);
        C_les = compute_convolution(uk_les, k_vec, K);

        % Dissipation terms
        D_dns = -(k_vec.^2) * nu .* uk_dns;
        nu_eff_les = nu + nu_t;
        D_les = -(k_vec.^2) .* nu_eff_les .* uk_les;

        % Explicit Euler integration
        uk_dns = uk_dns + dt * (-C_dns + D_dns);
        uk_les = uk_les + dt * (-C_les + D_les);

        % Apply large-scale forcing
        uk_dns(pk0 + 1) = uk1_const;
        uk_dns(pk0 - 1) = conj(uk1_const);
        uk_les(pk0 + 1) = uk1_const;
        uk_les(pk0 - 1) = conj(uk1_const);

        % Store total DNS energy
        energy_total_dns{idx}(t) = sum(abs(uk_dns).^2);

        % Progress display
        if usr.displayProgress
            if mod(t, round(nt/10)) == 0
                fprintf('N=%d --> %.0f%%\n', N, 100*t/nt);
            end
        end
    end

    % Compute final energy spectra (positive wavenumbers only)
    k_pos = k_vec(pk0+1:end);
    E_k_dns = abs(uk_dns(pk0+1:end)).^2;
    E_k_les = abs(uk_les(pk0+1:end)).^2;

    % Plot energy spectra
    if usr.doPlots
        loglog(k_pos, E_k_dns, 'o-', 'LineWidth', 1.5, 'DisplayName', ['DNS N=' num2str(N)]);
        loglog(k_pos, E_k_les, 's-', 'LineWidth', 1.5, 'DisplayName', ['LES N=' num2str(N)]);
    end
    % Store time vector for later plotting
    time_all{idx} = dt * (0:nt-1);
end

if usr.doPlots
    %% 3 Energy Spectrum Plot (Figure 1)
    % Add reference slope for comparison
    k_ref = k_pos(1);
    ref = (k_pos / k_ref).^(-2) * E_k_dns(1);
    loglog(k_pos, ref, '--k', 'LineWidth', 1.2, 'DisplayName', 'Slope -2 Reference');
    
    xlabel('Wavenumber k');
    ylabel('Energy E(k)');
    title('Energy Spectra Comparison (DNS vs LES)');
    legend('Location', 'southwest');
    hold off;
    
    %% 4 Total Energy Evolution (Figure 2)
    figure(2);
    hold on;
    for idx = 1:length(N_list)
        plot(time_all{idx}, energy_total_dns{idx}, 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_list(idx))]);
    end
    xlabel('Time t');
    ylabel('Total Energy');
    title('Total DNS Energy Evolution for Different N');
    grid on;
    legend('Location', 'best');
    hold off;
end