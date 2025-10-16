clear all;
close all;
clc;

% Author: Joel Campo
% 220311 - Computational Engineering. MUEA
% Assignment 2 (CFD): Generic convection-diffusion transport equation

%% INPUT DATA:
% Case study differentation & general inputs
case_study = 1; % 1 for Updwind; 2 for Exponential; 3 for CDS; 4 for HDS; 5 for PLDS
phi_inlet = 1; % Generic phi inlet
L = 10; % [m]

diag = 1;
if diag == 1
    case_study = 1;
    v_0 = 10;
end

% Fluid properties
rho = 1.225; % kg/m^3
GAMMA = [10 10^3 10^6];
gamma = rho./GAMMA; % revisar que ho fa bé

% Performance values
N = 500;
tolerance = 1e-06;
maxIt = 10^6;
time_step = 0.1;

%% PREVIOUS CALCULATIONS
H = L/2;
M = H / L * N;

if diag == 1
    H = L;
    M = N;
end

x_in = -L/2;
y_in = H;
x_out = L/2;
y_out = 0;

A_x = L/N;
A_y = H/M;

%% MESH GENERATION
center_x = (0:N) * A_x;
center_y = (0:M) * A_y;

nodvec_x = ((0:N-1) * A_x) + (A_x / 2);
nodvec_y = ((0:M-1) * A_y) + (A_y / 2);

nod_x = [0, nodvec_x, L];
nod_y = [0, nodvec_y, H];

% Preallocate
sx_w = zeros(size(nod_x));
sx_e = zeros(size(nod_x));
sx_p = zeros(size(nod_x));
sy_n = zeros(size(nod_y));
sy_s = zeros(size(nod_y));
sy_p = zeros(size(nod_y));

% X calculations
sx_w(2:end-1) = nod_x(2:end-1) - nod_x(1:end-2);
sx_e(2:end-1) = nod_x(3:end) - nod_x(2:end-1);

if numel(center_x) ~= numel(nod_x)
    center_x = linspace(min(nod_x), max(nod_x), numel(nod_x));
end
sx_p(2:end-1) = center_x(2:end-1) - center_x(1:end-2);

% X boundary conditions
sx_w(1) = 0;
sx_w(end) = nod_x(end) - nod_x(end-1);
sx_e(1) = nod_x(2) - nod_x(1);
sx_e(end) = 0;
sx_p(1) = center_x(2) / 2;
sx_p(end) = (center_x(end) - center_x(end-1)) / 2;

% Y calculations
sy_s(2:end-1) = nod_y(2:end-1) - nod_y(1:end-2);
sy_n(2:end-1) = nod_y(3:end) - nod_y(2:end-1);

if numel(center_y) ~= numel(nod_y)
    center_y = linspace(min(nod_y), max(nod_y), numel(nod_y));
end
sy_p(2:end-1) = center_y(2:end-1) - center_y(1:end-2);

% Y boundary conditions
sy_s(1) = 0;
sy_s(end) = nod_y(end) - nod_y(end-1);
sy_n(1) = nod_y(2) - nod_y(1);
sy_n(end) = 0;
sy_p(1) = center_y(2) / 2;
sy_p(end) = (center_y(end) - center_y(end-1)) / 2;

% Meshgrid preparation
x = nod_x - H;

[surf_x, surf_y] = meshgrid(x, fliplr(nod_y));
surf_x = surf_x';
surf_y = fliplr(surf_y');

%% Velocities computation
v_x = 2 * surf_y .* (1 - surf_x.^2);
v_y = -2 * surf_x .* (1 - surf_y.^2);

if diag == 1
    v_x(:,:) = v_0 * cosd(45);
    v_y(:,:) = v_0 * sind(45);
end

%% SOLVER
% Initialitation
phi_in = ones(N+2, M+2) * phi_inlet;
phi_prev = phi_in;
phi = zeros(N+2, M+2);
phi_def = zeros(N+2, M+2, length(gamma));

% Calculations
for l = 1:length(gamma) 
    a_n = zeros(N+2,M+2);
    a_s = zeros(N+2,M+2);
    a_e = zeros(N+2,M+2);
    a_w = zeros(N+2,M+2);
    a_p = ones(N+2,M+2);
    a_p_base = zeros(N+2,M+2);
    b_p = zeros(N+2,M+2);

    % Efficiency
    D_n = (gamma(l) * A_x) ./ sy_n(2:M+1);
    D_s = (gamma(l) * A_x) ./ sy_s(2:M+1);
    D_e = (gamma(l) * A_y) ./ sx_e(2:N+1);
    D_w = (gamma(l) * A_y) ./ sx_w(2:N+1);

    for j = 2:M+1
        for i = 2:N+1
            % Fluxes
            F_n = rho * (v_y(i,j+1) + v_y(i,j)) * A_x / 2;
            F_s = rho * (v_y(i,j-1) + v_y(i,j)) * A_x / 2;
            F_e = rho * (v_x(i+1,j) + v_x(i,j)) * A_y / 2;
            F_w = rho * (v_x(i-1,j) + v_x(i,j)) * A_y / 2;

            % Parameters
            p_n = F_n ./ D_n(j-1);
            p_s = F_s ./ D_s(j-1);
            p_e = F_e ./ D_e(i-1);
            p_w = F_w ./ D_w(i-1);

            if case_study == 1 % Upwind
                A_n = 1;
                A_s = 1;
                A_e = 1;
                A_w = 1;

            elseif case_study == 2 % Exponential
                A_n = abs(p_n) / exp(abs(p_n) - 1);
                A_s = abs(p_s) / exp(abs(p_s) - 1);
                A_e = abs(p_e) / exp(abs(p_e) - 1);
                A_w = abs(p_w) / exp(abs(p_w) - 1);

            elseif case_study == 3 % CDS
                A_n = abs((abs(p_n)/2)-1);
                A_s = abs((abs(p_s)/2)-1);
                A_e = abs((abs(p_e)/2)-1);
                A_w = abs((abs(p_w)/2)-1);

            elseif case_study == 4 % HDS
                A_n = max(0, 1 - (abs(p_n)/2));
                A_s = max(0, 1 - (abs(p_s)/2));
                A_e = max(0, 1 - (abs(p_e)/2));
                A_w = max(0, 1 - (abs(p_w)/2));

            elseif case_study == 5 % PLDS
                A_n = max(0, 1 - (abs(p_n)/2)^5);
                A_s = max(0, 1 - (abs(p_s)/2)^5);
                A_e = max(0, 1 - (abs(p_e)/2)^5);
                A_w = max(0, 1 - (abs(p_w)/2)^5);

            else
                error("The case_study selected is not available");
            end
 
            a_n(i,j) = D_n(j-1) * A_n + max(-F_n, 0);
            a_s(i,j) = D_s(j-1) * A_s + max(F_s, 0);
            a_e(i,j) = D_e(i-1) * A_e + max(-F_e, 0);
            a_w(i,j) = D_w(i-1) * A_w + max(F_w, 0);

            a_p_base(i,j) = rho * sx_p(i) * sy_p(j) / time_step;
            a_p(i,j) = a_e(i,j) + a_w(i,j) + a_n(i,j) + a_s(i,j) + a_p_base(i,j);
        end
    end
    % Boundaries
    b_p(1, 2:M+1) = 1 - tanh(10);
    b_p(N+2, 2:M+1) = 1 - tanh(10);
    b_p(2:N+1, M+2) = 1 - tanh(10);

    for i = 2:N+1
        if surf_x(i, M+2) < 0
            b_p(i, 1) = 1 + tanh(10 * (2 * surf_x(i, M+2) + 1));
        else
            b_p(i, 1) = a_n(i, 1) + a_s(i, 1) + a_e(i, 1) + a_w(i, 1) +...
                a_p_base(i,1);
            a_n(i, 1) = 1;
        end
    end

    if diag == 1
        b_p(1,1:M+2) = 100;
        b_p(N+2,1:M+2) = 0;
        b_p(1:N+2,1) = 0;
        b_p(1:N+2,M+2) = 100;
    end

    for k = 1:maxIt 
        b_p(2:N+1, 2:M+1) = a_p_base(2:N+1, 2:M+1) .* phi_prev(2:N+1, 2:M+1);
        coef.bp = b_p;

        for m = 1:maxIt
            phi = computePhi(N, M, a_n, a_s, a_e, a_w, a_p, b_p, phi_in);
            err_phi = max(abs(phi - phi_in), [], 'all');
            phi_in = phi;
            if err_phi < tolerance
                break
            end
        end

        err = abs(phi_prev - phi);
        error_value = max(max(err./time_step));
        phi_prev = phi;
        if error_value <= tolerance
            break
        end

    end

    phi_def(:,:,l) = phi; % Data storage
end

%% RESULTS VERIFICATION
if diag == 0
r = [0.000 0.100 0.200 0.300 0.400 0.500 0.600 0.700 0.800 0.900 1.000];
p = [1.989 1.402 1.146 0.946 0.775 0.621 0.480 0.349 0.227 0.111 0.000;
    2.00 1.999 1.9997 1.985 1.841 0.951 0.154 0.001 0.000 0.000 0.000;
    2.000 2.000 2.000 1.999 1.964 1.000 0.036 0.001 0.000 0.000 0.000];
r = transpose(r);
p = transpose(p);

phi_n = zeros(length(r),size(p,1));

for i = 1:size(p,2)
    assert(r(1) >= (nod_x(1)-H), 'You must have min(x) <= min(xx).');
    assert(r(end) <= (nod_x(end)-H), 'You must have max(xx) <= max(x).');
    
    int_phi = interp1(nod_x - H, phi_def(:,1,i), r, 'linear', 'extrap');
    phi_n(:,i)= int_phi(:);
end

figure;
hold on
for i = 1:length(gamma)
    plot(r, phi_n(:,i));
    plot(r, p(:,i), 'Marker', 'o', 'MarkerSize', 5);
end
grid on;
xlabel('x [m]');
ylabel('\rho/\Gamma');
ylim([-0.1, 2]);
title('Code verification (analytic vs numeric)', 'FontSize', 14);
legend('\rho/\Gamma = 10','Analytic \rho/\Gamma = 10', ...
    '\rho/\Gamma = 10^3','Analytic \rho/\Gamma = 10^3', ...
    '\rho/\Gamma = 10^6','Analytic \rho/\Gamma = 10^6');
end

%% PLOTTING
if diag == 0
    for i = 1:length(gamma)
        figure(i+1)
        surf(surf_x, surf_y, phi_def(:,:,i), 'Edgecolor', 'none', 'FaceColor', 'interp');
        c = colorbar;
        c.Label.String = '\Phi';
        colormap jet
        title('Distribution map \rho/\Gamma = ', num2str(gamma(i)))
        view(0, 90);
        axis equal
        xlabel('x [m]')
        ylabel('y [m]')
        xlim([-1, 1])
        ylim([0, 1])
    end
elseif diag == 1
    for i = 1:length(gamma)
        figure(i)
        surf(surf_x, surf_y, phi_def(:,:,i), 'Edgecolor', 'none', 'FaceColor', 'interp');
        view(0,90)
        c = colorbar;
        colormap jet
        grid off
        title('Distribution map \rho/\Gamma = ', num2str(gamma(i)))
        axis equal
        xlabel('x [m]')
        ylabel('y [m]')
        xlim([-10, 0])
        ylim([0, 10])
    end
end