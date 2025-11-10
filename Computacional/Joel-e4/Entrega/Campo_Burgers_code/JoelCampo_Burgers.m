clear all;
close all;
clc;

% Author: Joel Campo
% 220311 - Computational Engineering. MUEA
% Assignment 4 (CFD): Turbulence

%% INPUT DATA:
% Physical data
Re = 40; % Reynolds number
C_k = 0.4523; % Kolmogorov values

% Other data
C1 = 0.001;
N = 20; % Fourier modes
m = 2; % Slope in the energy plot

% Performance values
tolerance = 1e-05;
maxIt = 10^5;

%% PREVIOUS CALCULATIONS
k = 1:N; % Mode number
u = 1./k; % Velocities
time_step = C1 * Re/(N^2);

E_k_store = zeros(1,N,3);

%% SOLVER
for case_study = 1:3 % 1 for DNS, 2+ for LES with combinations
    for i = 1:maxIt
        % Preallocate
        conv = zeros(1,N);
        diff = zeros(1,N);

        for j = 1:N
            % Convective term
            for p = (-N+j):N
                q = j-p;
                if p~=0 && q~=0
                    if p > 0
                        u_p = u(p);
                    else
                        u_p = conj(u(-p));
                    end
    
                    if q > 0
                        u_q = u(q);
                    else
                        u_q = conj(u(-q));
                    end
                    conv(j) = conv(j) + u_p * q * 1i * u_q;
                end
            end

            % Diffusive term
            if case_study == 1 % DNS
                diff(j) = (k(j)^2) * u(j) / Re;

            elseif case_study == 2 % LES
                MetLes = 0.31 * (5 - m)/(m + 1) * sqrt(3 - m) * C_k^(-3/2);
                eddy = 1 + 34.5 * exp(-3.03 * (k(N)/k(j)));
                E_k_n = u(N) * conj(u(N));
                nu = MetLes * (E_k_n / k(N))^0.5 * eddy;

                diff(j) = k(j)^2 * u(j) * (1/Re + nu);

            elseif case_study == 3 % LES C_k = 0.05
                C_k = 0.05;
                MetLes = 0.31 * (5 - m)/(m + 1) * sqrt(3 - m) * C_k^(-3/2);
                eddy = 1 + 34.5 * exp(-3.03 * (k(N)/k(j)));
                E_k_n = u(N) * conj(u(N));
                nu = MetLes * (E_k_n / k(N))^0.5 * eddy;

                diff(j) = k(j)^2 * u(j) * (1/Re + nu);
            else
                error('The case study does not exist');
            end
        end

        % Velocities
        u_after = u + time_step * (conv - diff);
        u_after(1) = 1;

        % Error
        err = abs(u_after - u);
        if max(err(:)) < tolerance
            break
        else
            u = u_after;
        end
    end

    % Kinetic energy
    E_k = u_after.*conj(u_after);
    E_k_store(:,:,case_study) = E_k;
end

if Re == 20
    DNS_load = load("DNS20.mat");
elseif Re == 40
    DNS_load = load("DNS40.mat");
elseif Re == 80
    DNS_load = load("DNS80.mat");
elseif Re == 100
    DNS_load = load("DNS100.mat");
elseif Re == 180
    DNS_load = load("DNS180.mat");
else
    error("DNS(N=100) not computed for this Reynolds number");
end

figure(1)
loglog (k,E_k_store(:,:,3),'r','Marker','.')
hold on
loglog (k,E_k_store(:,:,2),'m','Marker','.')
hold on
loglog (k,E_k_store(:,:,1),'b','Marker','.')
hold on
k = 1:100;
loglog (k,DNS_load.E_k,'g','Marker','.')
hold on
loglog (k,1./k.^m, 'k-.', 'LineWidth', 2)
legend ('LES (N = 20; C_k = 0.4523)','LES (N = 20; C_k = 0.0500)', ...
    'DNS (N = 20)','DNS (N = 100)','Slope (m = -2)')
xlabel('k')
ylabel('Kinetic energy ( E_k )')
title('Energy spectrum for Re = 100')