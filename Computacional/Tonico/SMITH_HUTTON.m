% ANTONIO LUQUE LÓPEZ
% MUEA 
% ASSIGNMENT 2: SMITH-HUTTON RESOLUTION
%--------------------------------------------------------------------------

clear workspace,
clear all;
clc;        
clearvars;  
close all;  
%--------------------------------------------------------------------------

%%INSTRUCTIONS%%

%To plot Phi(x,y)-->plotPhi=1. If not, plotPhi=0(Recommended)
%For rg=[50, 75, 100]-->rgval=2. For rg=[10, 1e3, 1e6]-->rgval=1.(CDS scheme is not valid for high ρ/Γ-->Code crashes in line 68 becouse the "for" visits the scheme 2(CDS). To use rgval=1, change the "1:5" for a value different to 2.)

PlotPhi=0; rgval=2;

%--------------------------------------------------------------------------

% 1. INPUT DATA
% 1.1 Physical data
L=2;
H=1;
alpha=10;

% 1.2 Numerical data
N = 100;       
M = 50; 
delta=1e-10;
scheme=1;

% 2. PREVIOUS CALCULATIONS
dx = L/N;   
dy = H/M;  

[x, y] = SH_mesh_generator(N,M,L,H); % Function that creates the geometry of the mesh.
[X, Y] = meshgrid(x, y);

R=287;
P_ref=101325;
T_ref=300;

rho=ones(M+2,N+2)*(P_ref/(R*T_ref));

if rgval==1
rg = [10, 1e3, 1e6]; % ρ/Γ values 1
else
rg = [50, 75, 100]; % ρ/Γ values 2
end

CTTC_data = [ % values given as a solution for Phi(ρ/Γ=[10, 1e3, 1e6])
    1.9890, 2.0000, 2.0000;
    1.4020, 1.9990, 2.0000;
    1.1460, 1.9997, 2.0000;
    0.9460, 1.9850, 1.9990;
    0.7750, 1.8410, 1.9640;
    0.6210, 0.9510, 1.0000;
    0.4800, 0.1540, 0.0360;
    0.3490, 0.0010, 0.0010;
    0.2270, 0.0000, 0.0000;
    0.1110, 0.0000, 0.0000;
    0.0000, 0.0000, 0.0000];

x_output = 0:0.1:1;
phi_values = zeros(length(x_output), length(rg));  % almacenar interpolaciones
phi_scheme=zeros(length(x_output),5*length(rg));

% 3. COMPUTE PHI MATRIX 
for scheme=1:5 %Loop for plotting the differents phi values depending on the numerical scheme for each rho/Gamma values
    for k = 1:length(rg)%Loop for creating the Phi(x) table depending on rho/Gamma values
        rgi = rg(k);
        Gamma = (rho / rgi) .* ones(M+2, N+2);
    
        [u,v] = set_velocity_field(N, M, x, y); % Creates the velocity field
        [a_E, a_W, a_N, a_S, a_P, b, AP_e, AP_w, AP_n, AP_s] = compute_coefficients(u, v, rho, Gamma, x, y, dx, dy, N, M, scheme); % Compute the discretization coefficients
        [Phi, iter, error] = SH_gauss_seidel(a_E, a_W, a_N, a_S, a_P, b, delta, M, N, alpha, x); % Solve the discretized equation
    
        % Creates the Phi interpolation for x=[0, 0.1, 0.2, ..., 0.9, 1]
        a = ceil((N+2)/2+1) : N+2;
        x_cols = x(a);
        Phi_output = Phi(1, a);
        Phi_interp = interp1(x_cols, Phi_output, x_output, 'makima');
    
        phi_values(:, k) = Phi_interp;

            %phi_scheme creates a matrix with Phi for each scheme for each rho/gamma
            c = (scheme-1)*length(rg) + k;
            phi_scheme(:, c) = Phi_interp;

        %Plot Phi(x,y)
        if PlotPhi==1
        figure;
        pcolor(X, Y, Phi);
        shading interp;
        colorbar;
        title(sprintf('\\phi distribution – Scheme %d, \\rho/\\Gamma = %.0e', scheme, rgi));
        xlabel('x'); ylabel('y');
        end
    end  
end

table_rg(x_output, phi_values, rg); %Function that prints the table with Phi for each rho/gamma


err=abs(phi_values-CTTC_data); %Compute the absolute error of Phi

%Plot a graphic for comparing eache scheme for each rho/gamma
if rgval==2 
figure;
plot(x_output, phi_scheme(:, [1, 4, 7, 10, 13]), 'LineWidth', 1.5)
title('Comparison of schemes for \rho/\Gamma = 50')
xlabel('x-position')
ylabel('\phi at y = 0')
legend('UDS', 'CDS', 'HDS', 'EDS', 'PLDS', 'Location', 'best');

figure;
plot(x_output, phi_scheme(:, [2, 5, 8, 11, 14]), 'LineWidth', 1.5)
title('Comparison of schemes for \rho/\Gamma = 75')
xlabel('x-position')
ylabel('\phi at y = 0')
legend('UDS', 'CDS', 'HDS', 'EDS', 'PLDS', 'Location', 'best');

figure;
plot(x_output, phi_scheme(:, [3, 6, 9, 12, 15]), 'LineWidth', 1.5)
title('Comparison of schemes for \rho/\Gamma = 100')
xlabel('x-position')
ylabel('\phi at y = 0')
legend('UDS', 'CDS', 'HDS', 'EDS', 'PLDS', 'Location', 'best');

end
