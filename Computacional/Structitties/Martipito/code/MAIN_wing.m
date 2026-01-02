%% PART III - WING MODELLING (TEMPLATE)
% Useful commands to initialize script
clear
close all
clc

%% SOLUTION CONTROLS
problem_type = 3;  % 1: Bending, 2: Torsion, 3: Pressure + Weight
solution_type = 1; % 1: Direct Method, 2: Model Order Reduction (MOR)
Nm = 10;           % Number of modes to compute
Im = 1:9;          % Modes to use for MOR projection

% Component flags (1: Include, 0: Neglect)
WB = 1; % Wingbox
ST = 1; % Stringers
RB = 1; % Ribs
SK = 1; % Skin

%% DATA
% Material Properties
E = 75e9;       % Pa
v_poisson = 0.3;
rho = 2200;     % kg/m3
G = E/(2*(1+v_poisson));

% Geometry
c = 0.85; %m
b = 3; %m
y0 = 0.924; %m
y1 = 0.125; %m
y2 = 0.475; %m
yc = 0.285; %m

% Loads Data
F = 1; % N (Unit force for bending)
T = 1; % Nm (Unit torque)
p_inf = 12500;   % Pa (Pressure)
alpha = 10*pi/180; % rad (Angle of attack)
g0 = 9.80665;    % m/s^2

%% PREPROCESS
% Load mesh data
if isfile('wing.mat')
    load('wing.mat','xn','Tn_st','Tm_st','Tn_wb','Tm_wb','Tn_rb','Tm_rb','Tn_sk','Tm_sk','indRoot','indPointA','indPointB','indSpar1','indSpar2','n_u','n_l');
else
    error('File wing.mat not found.');
end

N_nodes = size(xn,1);
Ndof = 6*N_nodes;

% --- 1. Boundary Conditions (Up) ---
% Fix all DOFs at root nodes (Vectorized)
num_root_nodes = length(indRoot);
Up = zeros(num_root_nodes * 6, 3);
Up(:, 1) = 0;                                  % Displacement value
Up(:, 2) = repelem(indRoot, 6);                % Node indices
Up(:, 3) = repmat((1:6)', num_root_nodes, 1);  % DOFs 1 to 6

% --- 2. Load Definitions (Fe, Pe, Be) ---
% Calculate Load magnitudes
Fa_shear = F*((y2-yc)/(y2-y1));
Fb_shear = F*((yc-y1)/(y2-y1));
Fa_torque = -T/(y2-y1);
Fb_torque = T/(y2-y1);

% Initialize empty arrays
Fe_wb = []; Pe_wb = []; 
Fe_st = []; Qe_st = [];
Fe_rb = []; Pe_rb = [];
Fe_sk = []; Pe_sk = [];
Be = [];

% Define specific loads based on problem type
switch(problem_type)
    case 1 % Bending (Shear Force)
        Fe_wb = [Fa_shear indPointA 3;
                 Fb_shear indPointB 3];
    case 2 % Torsion (Torque)
        Fe_wb = [Fa_torque indPointA 3;
                 Fb_torque indPointB 3];
    case 3 % Pressure + Gravity
        % Pressure on Skin
        Pe_sk = ComputePressureState(xn,n_u,n_l,b,c,p_inf,alpha);
        % Body Force (Gravity) on all nodes
        node_indices = (1:N_nodes)';
        Be = [-g0 * ones(N_nodes,1), node_indices, 3*ones(N_nodes,1)];
end

% --- 3. Element Properties ---
% Wing Box (Shells)
h_vector_wb = [17.5e-3; 15e-3; 6e-3]; 
E_vector_wb = repmat(E,3,1);
v_vector_wb = repmat(v_poisson,3,1);
rho_vector_wb = repmat(rho,3,1);

% Stringers (Beams)
r_st = 5e-3; 
A_st = pi*r_st^2;
I_st = (pi/4)*r_st^4;
J_st = (pi/2)*r_st^4;
% Beam properties vectors
E_vector_st = E; G_vector_st = G; rho_vector_st = rho;
A_vector_st = A_st; Iy_vector_st = I_st; Iz_vector_st = I_st; J_vector_st = J_st;
ky_vector_st = 5/6; kz_vector_st = 5/6; kt_vector_st = 1;
j_orientations = [0 1 0]; 

% Ribs (Shells)
h_vector_rb = 5e-3; 
E_vector_rb = E; v_vector_rb = v_poisson; rho_vector_rb = rho;

% Skin (Shells)
h_vector_sk = 2e-3;
E_vector_sk = E; v_vector_sk = v_poisson; rho_vector_sk = rho;

%% SOLVER
% Obtain system matrices
if 1 % Set to 1 to (re)compute, 0 to load
    fprintf('Computing Wing Box Matrices...\n');
    [K_wb,M_wb,f_wb,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb,R_wb] = ComputeKM_shells(xn,N_nodes,size(Tn_wb,1),Tn_wb,Tm_wb,h_vector_wb,E_vector_wb,v_vector_wb,rho_vector_wb,Fe_wb,Be,Pe_wb);
    
    fprintf('Computing Stringer Matrices...\n');
    [K_st,M_st,f_st] = ComputeKM_beams(xn,N_nodes,size(Tn_st,1),Tn_st,Tm_st,E_vector_st,G_vector_st,rho_vector_st,j_orientations,A_vector_st,Iy_vector_st,Iz_vector_st,J_vector_st,ky_vector_st,kz_vector_st,kt_vector_st,Fe_st,Be,Qe_st);
    
    fprintf('Computing Rib Matrices...\n');
    [K_rb,M_rb,f_rb,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb,R_rb] = ComputeKM_shells(xn,N_nodes,size(Tn_rb,1),Tn_rb,Tm_rb,h_vector_rb,E_vector_rb,v_vector_rb,rho_vector_rb,Fe_rb,Be,Pe_rb);
    
    fprintf('Computing Skin Matrices...\n');
    [K_sk,M_sk,f_sk,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk,R_sk] = ComputeKM_shells(xn,N_nodes,size(Tn_sk,1),Tn_sk,Tm_sk,h_vector_sk,E_vector_sk,v_vector_sk,rho_vector_sk,Fe_sk,Be,Pe_sk);
    
    % Save results
    save('wing_matrices.mat','K_wb','M_wb','K_st','M_st','K_rb','M_rb','K_sk','M_sk','f_wb','f_st','f_rb','f_sk');
    save('wing_aux_matrices.mat','Bb_wb','Bmn_wb','Bmt_wb','Bs_wb','R_wb','Bb_rb','Bmn_rb','Bmt_rb','Bs_rb','R_rb','Bb_sk','Bmn_sk','Bmt_sk','Bs_sk','R_sk');
    fprintf('Matrices saved.\n');
else
    load('wing_matrices.mat');
    load('wing_aux_matrices.mat');
    fprintf('Matrices loaded.\n');
end

% Assemble Global System
fprintf('Assembling Global System...\n');
K = WB*K_wb + ST*K_st + RB*K_rb + SK*K_sk;
M = WB*M_wb + ST*M_st + RB*M_rb + SK*M_sk;
f = f_wb + f_st + f_rb + f_sk;

% Partition System (Free vs Prescribed DOFs)
prescribed_dofs = 6*(Up(:,2)-1) + Up(:,3);
free_dofs = setdiff(1:Ndof, prescribed_dofs);

% --- Modal Analysis ---
fprintf('Performing Modal Analysis...\n');
% Using MATLAB eigs instead of custom MAIN_modal
[V, D] = eigs(K(free_dofs,free_dofs), M(free_dofs,free_dofs), Nm, 'sm');
[lambda_sorted, idx] = sort(diag(D));
V_sorted = V(:, idx);
freq = sqrt(lambda_sorted)/(2*pi);

Phi = zeros(Ndof, Nm);
Phi(free_dofs, :) = V_sorted;

% Normalize Mass (Optional but recommended)
for k = 1:Nm
    norm_factor = sqrt(Phi(:,k)' * M * Phi(:,k));
    Phi(:,k) = Phi(:,k) / norm_factor;
end
fprintf('Natural frequencies computed. First freq: %.2f Hz\n', freq(1));

% --- Static Solution ---
fprintf('Solving Static System...\n');
u = zeros(Ndof, 1);
u(prescribed_dofs) = Up(:,1);

if solution_type == 1 % Direct Method
    % Ku = f -> u_f = K_ff \ (f_f - K_fp * u_p)
    f_eff = f(free_dofs) - K(free_dofs, prescribed_dofs) * u(prescribed_dofs);
    u(free_dofs) = K(free_dofs, free_dofs) \ f_eff;
    
elseif solution_type == 2 % Model Order Reduction
    % u ~ Sum(Phi_i * alpha_i)
    % alpha_i = (Phi_i' * F) / (lambda_i - omega^2) -> For static, omega=0
    
    % We need enough modes for accuracy. Check if Nm is sufficient.
    u_reduced = zeros(Ndof, 1);
    for k = Im
        modal_force = Phi(:,k)' * f;
        stiffness_val = lambda_sorted(k); % eigenvalue = omega_n^2
        alpha = modal_force / stiffness_val;
        u_reduced = u_reduced + Phi(:,k) * alpha;
    end
    u = u_reduced; % Ignoring boundary corrections for simplicity in MOR
end

% --- Stress Calculation ---
fprintf('Calculating Stresses...\n');
% We calculate stresses for shell components
sigVM_wb = Compute_stress_shells(size(Tn_wb,1),u,Tn_wb,Tm_wb,h_vector_wb,v_vector_wb,E_vector_wb,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb,R_wb);
sigVM_rb = Compute_stress_shells(size(Tn_rb,1),u,Tn_rb,Tm_rb,h_vector_rb,v_vector_rb,E_vector_rb,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb,R_rb);
sigVM_sk = Compute_stress_shells(size(Tn_sk,1),u,Tn_sk,Tm_sk,h_vector_sk,v_vector_sk,E_vector_sk,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk,R_sk);

%% POSTPROCESS
save('wing_results.mat');

% --- 1. 2D Plots: Spanwise Deflection & Twist ---
n_points = size(indSpar1,1);
theta_x = zeros(n_points,1);
u_z = zeros(n_points,1);
xn_plot = xn(indSpar1,1);

for i = 1:n_points
    node_front = indSpar1(i);
    node_rear = indSpar2(i);
    
    uz_front = u(6*(node_front-1) + 3);
    uz_rear  = u(6*(node_rear-1) + 3);
    
    % Twist = (delta_z) / distance
    dist_y = y2 - y1; % Approximate, or calculate norm(xn(node_rear)-xn(node_front))
    theta_x(i) = (uz_rear - uz_front) / dist_y;
    
    % Uz at Shear Center
    u_z(i) = uz_front + theta_x(i)*(yc - y1);
end

figure('Name', 'Spanwise Results', 'Color', 'w');

subplot(2,1,1);
plot(xn_plot, u_z, 'LineWidth', 1.5);
title('Vertical Deflection (u_z) along Span');
xlabel('Spanwise Coordinate x [m]');   % <--- Unidades añadidas
ylabel('Vertical Deflection u_z [m]'); % <--- Unidades añadidas
grid on;
legend(sprintf('Problem Type %d', problem_type));

subplot(2,1,2);
plot(xn_plot, theta_x, 'LineWidth', 1.5);
title('Twist Angle (\theta_x) along Span');
xlabel('Spanwise Coordinate x [m]'); % <--- Unidades añadidas
ylabel('Twist Angle \theta_x [rad]'); % <--- Unidades añadidas
grid on;

% --- 2. 3D Deformation & Stress Plots ---
% Determine scale factor automatically if needed
scale = 0.1 * max(abs(xn(:,1))) / max(abs(u(3:6:end)) + eps);
fprintf('Plotting deformations with scale factor: %.2f\n', scale);

% Se añaden unidades [Pa] a los títulos de los plots 3D
plotDeformed('wing', xn, Tn_wb, u, scale, sigVM_wb); 
title(['Wingbox Deformation & Von Mises Stress [Pa] - Scale: ' num2str(scale)]); % <--- Unidad añadida

plotDeformed('wing', xn, Tn_rb, u, scale, sigVM_rb); 
title(['Ribs Deformation & Von Mises Stress [Pa] - Scale: ' num2str(scale)]);    % <--- Unidad añadida

plotDeformed('wing', xn, Tn_sk, u, scale, sigVM_sk); 
title(['Skin Deformation & Von Mises Stress [Pa] - Scale: ' num2str(scale)]);    % <--- Unidad añadida

% --- 3. Mode Shapes ---
fprintf('Plotting Mode Shapes...\n');
plotModes('wing', Phi, freq, 1:min(9, Nm));

fprintf('--- Analysis Complete ---\n');