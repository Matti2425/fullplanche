%% PART I - BEAM MODELLING 
% Gerard Martínez Vilella
% David Rodríguez López

% Useful commands to initialize script
clear
close all
clearvars

%% DATA

problem_type = "torsion"; %1 for bending. 2 for torsion.
Nm = 6;

% Define the problem's data (e.g. dimensions, material properties, etc.)
E = 71e9; %Pa
v_poisson = 0.33;
G = E/(2*(1+v_poisson));
rho = 2300; %kg/m3
c = 0.8; %m
b = 5.5; %m
y0 = 0.3083; %m
y1 = 0.176; %m
y2 = 0.464; %m
h1 = 17.5e-3; %m
h2 = 15e-3; %m
h3 = 6e-3; %m
A = 0.0062; %m2
Iy = 8.801e-6; %m4
Iz = 8.069e-5; %m4
J = 8.949e-5; %m4
ky = 0.5457;
kz = 0.1767;
kt = 0.2920;
yc = 0.2954; %m

%% PREPROCESS
zeros(0,3);
% Load mesh data
load('beam.mat','xn','Tn','Tm');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 2]
% Tm : Material connectivities [Nelem x 1]

% Define boundary conditions and forces data matrices: Up, Fe, Qe, Be
N_elements = size(Tn,1);
N_nodes = size(xn,1);
Up = [0 1 1; % displacement, node, dof
      0 1 2;
      0 1 3;
      0 1 4;
      0 1 5;
      0 1 6];

% External forces (force, node, dof)
switch(problem_type)
    case "bending"
    Fe = [1 N_nodes 3];
    case "torsion"
    Fe = [1 N_nodes 5];
end

% Distributed forces (force, node, dof)
Qe = [0 1 2];

% Body forces (force, node, dof)
Be = [0 1 2];

j_orientations = [0 1 0]; 

%more rows would be added if more values existed
E_vector = [E]; 
A_vector = [A];  
G_vector = [G];  
Iy_vector = [Iy];  
Iz_vector = [Iz];  
ky_vector = [ky];  
kz_vector = [kz];  
kt_vector = [kt];  
J_vector = [J];  
rho_vector = [rho];  

%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
% ========================================================
%  SECTION 2 — MATRIX ASSEMBLY (sparse + blkdiag version)
% ========================================================

    Ndof = 6*N_nodes;
    K = sparse(Ndof,Ndof);
    M  = sparse(Ndof,Ndof);
    
    for e = 1:N_elements
    
        % --- Longitud y ejes locales ----------------------
        x1 = xn(Tn(e,1),:)';
        x2 = xn(Tn(e,2),:)';
        Le = norm(x2 - x1);
    
        i = (x2 - x1)/Le;
        j = j_orientations(Tm(e,1),:)';
        k = cross(i,j);
    
        % --- Rotación local Rp (6x6) ----------------------
        R3 = [i j k];
        Rp = blkdiag(R3, R3);
    
        % --- Rotación global (12x12) ----------------------
        Re = blkdiag(Rp, Rp);
    
        % --- Axial ----------------------------------------
        Nx = [-1/Le; 1/Le];
        Ba_row = [Nx(1) 0 0 0 0 0 Nx(2) 0 0 0 0 0];
        Ba = sparse(Ba_row);
    
        Ca = E_vector(Tm(e,1))*A_vector(Tm(e,1));
    
        Ka = Le * (Re' * (Ba' * Ca * Ba) * Re);
    
        % --- Flexión --------------------------------------
        Bb = sparse([ ...
            0 0 0 0 Nx(1) 0   0 0 0 0 Nx(2) 0;
            0 0 0 0 0 Nx(1)  0 0 0 0 0 Nx(2)]);
    
        Cb = E_vector(Tm(e,1)) * diag([Iy_vector(Tm(e,1)) Iz_vector(Tm(e,1))]);
    
        Kb = Le * (Re' * (Bb' * Cb * Bb) * Re);
    
        % --- Cortante --------------------------------------
        N_G = 1/2;
    
        Bs = sparse([ ...
            0 Nx(1) 0 0 0 -N_G  0 Nx(2) 0 0 0 -N_G;
            0 0 Nx(1) 0 N_G 0   0 0 Nx(2) 0 N_G 0 ]);
    
        Cs = G_vector(Tm(e,1))*A_vector(Tm(e,1))* ...
             diag([ky_vector(Tm(e,1)) kz_vector(Tm(e,1))]);
    
        Ks = Le * (Re' * (Bs' * Cs * Bs) * Re);
    
        % --- Torsión --------------------------------------
        Bt_row = [0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0 0];
        Bt = sparse(Bt_row);
    
        Ct = G_vector(Tm(e,1))*J_vector(Tm(e,1))*kt_vector(Tm(e,1));
    
        Kt = Le * (Re' * (Bt' * Ct * Bt) * Re);
    
        % -----------------------------------------------
        %     MASA CONSISTENTE (también local)
        % -----------------------------------------------
        Xi = [-1/sqrt(3); 1/sqrt(3)];
        w  = [1; 1];
    
        A  = A_vector(Tm(e,1));
        Iy = Iy_vector(Tm(e,1));
        Iz = Iz_vector(Tm(e,1));
        J  = J_vector(Tm(e,1));
        rho = rho_vector(Tm(e,1));
    
        rho_mat = rho * blkdiag(A*speye(3), J, Iy, Iz);  % 6×6
    
        Me = sparse(12,12);
    
        for kk = 1:2
            N = [(1 - Xi(kk))/2; (1 + Xi(kk))/2];
            Nrow = [N(1) 0 0 0 0 0 N(2) 0 0 0 0 0];
    
            Nmat = sparse(6,12);
                for k = 1:6
                    Nmat(k, k)     = N(1);
                    Nmat(k, 6 + k) = N(2);
                end

            Me = Me + (w(kk)*Le/2) * (Re' * (Nmat' * rho_mat * Nmat) * Re);
        end
    
        % --- DOFs globales --------------------------------
        Idof = zeros(12,1);
        for j = 1:6
            Idof(j)     = 6*(Tn(e,1)-1) + j;
            Idof(6 + j) = 6*(Tn(e,2)-1) + j;
        end
    
        % --- Ensamblaje global -----------------------------
        K(Idof,Idof) = K(Idof,Idof) + Ka + Kb + Ks + Kt;
        M(Idof,Idof)  = M(Idof,Idof)  + Me;
    
    end


    fprintf('--- Applying boundary conditions and loads ---\n');
    
    % --- External Force Vector ---
    F = zeros(Ndof, 1);
    F_node = Fe(1,2);
    F_dof = Fe(1,3);
    F_global = (F_node - 1) * 6 + F_dof;
    F(F_global) = Fe(1,1); % apply unit or specified load
    
    % --- Apply Boundary Conditions ---
    prescribed_dofs = (Up(:,2) - 1) * 6 + Up(:,3);
    free_dofs = setdiff(1:Ndof, prescribed_dofs);
    
    % --- Solve system ---
    U = zeros(Ndof, 1);
    U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);
    
    fprintf('--- Static solution complete ---\n');
    
    % Perform modal analysis
    fprintf('--- Performing modal analysis ---\n');
    n_modes = min(Nm, length(free_dofs));
    [Phi, Omega2] = eigs(K(free_dofs, free_dofs), M(free_dofs, free_dofs), n_modes, 'sm');
    frequencies = sqrt(diag(Omega2)) / (2*pi);
    
    fprintf('Natural frequencies [Hz]:\n');
    disp(frequencies');
    
end
%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('beam_results.mat');

% --- 1. Extract Displacements from Global Vector U ---
% U is stored as [ux1; uy1; uz1; thx1; thy1; thz1; ux2; ... ]
u_x = U(1:6:end);       % Axial displacement
u_y = U(2:6:end);       % Lateral displacement (Y)
u_z = U(3:6:end);       % Vertical displacement (Z)
theta_x = U(4:6:end);   % Torsion/Twist
theta_y = U(5:6:end);   % Bending rotation about Y
theta_z = U(6:6:end);   % Bending rotation about Z

% Define x-axis for plotting (Spanwise coordinate)
x_span = xn(:,1); 

% --- 2. Plot Translational Displacements (ux, uy, uz) ---
figure('Name', 'Translational Displacements', 'NumberTitle', 'off');

subplot(3,1,1);
plot(x_span, u_x, 'LineWidth', 1.5);
title('Axial Displacement (u_x)');
ylabel('u_x (m)'); grid on;

subplot(3,1,2);
plot(x_span, u_y, 'LineWidth', 1.5);
title('Lateral Displacement (u_y)');
ylabel('u_y (m)'); grid on;

subplot(3,1,3);
plot(x_span, u_z, 'LineWidth', 1.5);
title('Vertical Deflection (u_z)');
xlabel('Span x (m)'); ylabel('u_z (m)'); grid on;

% Add legend based on problem type
switch(problem_type)
    case "bending"
        sgtitle('Static Solution: Bending Test (Translations)');
    case "torsion"
        sgtitle('Static Solution: Torsion Test (Translations)');
end

% --- 3. Plot Rotational Displacements (theta_x, theta_y, theta_z) ---
figure('Name', 'Rotational Displacements', 'NumberTitle', 'off');

subplot(3,1,1);
plot(x_span, theta_x, 'LineWidth', 1.5);
title('Twist Angle (\theta_x)');
ylabel('\theta_x (rad)'); grid on;

subplot(3,1,2);
plot(x_span, theta_y, 'LineWidth', 1.5);
title('Rotation about Y (\theta_y)');
ylabel('\theta_y (rad)'); grid on;

subplot(3,1,3);
plot(x_span, theta_z, 'LineWidth', 1.5);
title('Rotation about Z (\theta_z)');
xlabel('Span x (m)'); ylabel('\theta_z (rad)'); grid on;

switch(problem_type)
    case "bending"
        sgtitle('Static Solution: Bending Test (Rotations)');
    case "torsion"
        sgtitle('Static Solution: Torsion Test (Rotations)');
end