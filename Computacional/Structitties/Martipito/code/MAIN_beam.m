%% PART I - BEAM MODELLING 
% Pau Cornudella Quer
% Martí Güell Domingo

clear
close all
clearvars

%% DATA

problem_type = "torsion"; % CAS 1: "bending". CAS 2: "torsion".
Nm = 6; % N vibration modes in the modal analysis

% Define the problem's data (e.g. dimensions, material properties, etc.)
E = 75e9; %Pa
v_poisson = 0.3;
G = E/(2*(1+v_poisson));
rho = 2200; %kg/m3
c = 0.85; %m
b = 3; %m
y0 = 0.294; %m
y1 = 0.125; %m
y2 = 0.475; %m
h1 = 35e-3; %m
h2 = 30e-3; %m
h3 = 25e-3; %m
A = 0.0246; %m2
Iy = 0.727e-4; %m4
Iz = 4.002e-4; %m4
J = 4.729e-4; %m4
ky = 0.6873;
kz = 0.1843;
kt = 0.4602;
yc = 0.285; %m

%% PREPROCESS
zeros(0,3);
% Load mesh data
load('beam.mat','xn','Tn','Tm');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 2] %Conection between Tn(e,1) to Tn(e,2)
% Tm : Material connectivities [Nelem x 1]

% Define boundary conditions and forces data matrices: Up, Fe, Qe, Be
N_elements = size(Tn,1); 
N_nodes = size(xn,1);

Up = [0 1 1; % displacement, node, dof
      0 1 2; 
      0 1 3;
      0 1 4;
      0 1 5;
      0 1 6]; % All 6 DOFs of node 1 = 0. Clamped

% External forces (force, node, dof)
switch(problem_type)
    case "bending"
    Fe = [1 N_nodes 3]; %Force = 1 at last node in DOF 3 (uz)
    case "torsion"
    Fe = [1 N_nodes 4]; %Force = 1 at last node in DOF 4 (theta_x)
end

% Distributed forces (force, node, dof)
Qe = [0 1 2];

% Body forces (force, node, dof)
Be = [0 1 2];

j_orientations = [0 1 0]; %input data, volem que sigui cap a y global

%For convenience;
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

    %Matrix assembly (sparse + blkdiag version)
    Ndof = 6*N_nodes; 
    K = sparse(Ndof,Ndof); % Stiffness matrix
    M  = sparse(Ndof,Ndof);% Mass matrix
    
    for e = 1:N_elements
    
        % Lenght and local axes
        x1 = xn(Tn(e,1),:)'; % Node 1 of element e
        x2 = xn(Tn(e,2),:)'; % Node 2 of element e
        Le = norm(x2 - x1);  % Element lenght
    
        i = (x2 - x1)/Le; % Local X axis unit vector
        j = j_orientations(Tm(e,1),:)'; % Local Y axis unit vector
        k = cross(i,j); %Local Z axis unit vector
    
        % Rotación local Rp (6x6) 
        R3 = [i j k];
        Rp = blkdiag(R3, R3);
    
        % Rotación global (12x12)
        Re = blkdiag(Rp, Rp);
    
        % Axial (pag 6)
        Nx = [-1/Le; 1/Le]; 
        Ba_row = [Nx(1) 0 0 0 0 0 Nx(2) 0 0 0 0 0]; % Al node
        Ba = sparse(Ba_row);
        Ca = E_vector(Tm(e,1))*A_vector(Tm(e,1));
        Ka = Le*(Re'*(Ba'*Ca*Ba)*Re);
        %K glb = ....(K local).......
    
        % Flexió (pag 6)
        Bb = sparse([0 0 0 0 Nx(1) 0   0 0 0 0 Nx(2) 0;
                     0 0 0 0 0 Nx(1)  0 0 0 0 0 Nx(2)]);
        Cb = E_vector(Tm(e,1)) * diag([Iy_vector(Tm(e,1)) Iz_vector(Tm(e,1))]);
        Kb = Le*(Re'*(Bb'*Cb*Bb)*Re);
    
        % Tallant (pag 7)
        N_G = 1/2;
        Bs = sparse([0 Nx(1) 0 0 0 -N_G  0 Nx(2) 0 0 0 -N_G;
                     0 0 Nx(1) 0 N_G 0   0 0 Nx(2) 0 N_G 0 ]);
        Cs = G_vector(Tm(e,1))*A_vector(Tm(e,1))*diag([ky_vector(Tm(e,1)) kz_vector(Tm(e,1))]);
        Ks = Le*(Re'*(Bs'*Cs*Bs)*Re);
    
        %Torsió (pag 7)
        Bt_row = [0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0 0];
        Bt = sparse(Bt_row);
        Ct = G_vector(Tm(e,1))*J_vector(Tm(e,1))*kt_vector(Tm(e,1));
        Kt = Le*(Re'*(Bt'*Ct*Bt)*Re);
    
        %Ass Martix (local) (pag 8)
        Xi = [-1/sqrt(3); 1/sqrt(3)];
        w  = [1; 1];
        
        A  = A_vector(Tm(e,1));
        Iy = Iy_vector(Tm(e,1));
        Iz = Iz_vector(Tm(e,1));
        J  = J_vector(Tm(e,1));
        rho = rho_vector(Tm(e,1));
    
        rho_mat = rho * blkdiag(A*speye(3), J, Iy, Iz);  % 6×6 speye=identitat (mantenim sparse)
    
        Me = sparse(12,12);
    
        for kk = 1:2
            N = [(1 - Xi(kk))/2; (1 + Xi(kk))/2];
            Nrow = [N(1) 0 0 0 0 0 N(2) 0 0 0 0 0];
            Nmat = sparse(6,12);
                for k = 1:6
                    Nmat(k, k)   = N(1);
                    Nmat(k, 6+k) = N(2);
                end

            Me = Me + (w(kk)*Le/2) * (Re'*(Nmat'*rho_mat*Nmat) * Re);
        end
    
        % DOFs globals (pag 9)
        Idof = zeros(12,1);
        for j = 1:6
            Idof(j) = 6*(Tn(e,1)-1) + j;
            Idof(6 + j) = 6*(Tn(e,2)-1) + j;
        end
    
        % Ensamblatge (pag 9)
        K(Idof,Idof) = K(Idof,Idof) + Ka + Kb + Ks + Kt;
        M(Idof,Idof)  = M(Idof,Idof)  + Me;
    
    end


    fprintf('--- Applying boundary conditions and loads ---\n');
    
    % External Force Vector (pag 10)
    F = zeros(Ndof, 1);
    F_node = Fe(1,2);
    F_dof = Fe(1,3);
    F_global = (F_node - 1) * 6 + F_dof; %nodee to global dofs
   
    F(F_global) = Fe(1,1); % Specified load
    
    % Apply Boundary Conditions (pag 12)
    prescribed_dofs = (Up(:,2) - 1) * 6 + Up(:,3);
    free_dofs = setdiff(1:Ndof, prescribed_dofs);
    
    % Solve system (pag 12)
    U = zeros(Ndof, 1);
    U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);
    
    fprintf('--- Static solution complete ---\n');
    
    % Perform modal analysis (pag 4 frequency)
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
u_y = U(2:6:end);       % Lateral displacement 
u_z = U(3:6:end);       % Vertical displacement 
theta_x = U(4:6:end);   % Torsion
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

% --- 4. Plot mode shapes
 
Phi_full = zeros(Ndof, n_modes);
Phi_full(free_dofs, :) = Phi;
x_span = xn(:,1);
%plotModes('Beam Vibration Modes', xn, Tn, Phi_full, frequencies, 1:n_modes);


    
for i = 1:n_modes
        
    uy = Phi_full(2:6:end, i);
    uz = Phi_full(3:6:end, i); 
    thx = Phi_full(4:6:end, i);

    % Normalize for visualization
    scale = max(abs([uy; uz; thx]));      
    uy = uy/scale;   
    uz = uz/scale;
    thx = thx/scale;

    figure;
    subplot(3,1,1)
    plot(x_span, uy, 'LineWidth', 1.5); grid on;
    ylabel('$u_y$', 'Interpreter','latex');

   subplot(3,1,2)
   plot(x_span, uz, 'LineWidth', 1.5); grid on;
   ylabel('$u_z$', 'Interpreter','latex');

   subplot(3,1,3)
   plot(x_span, thx, 'LineWidth', 1.5); grid on;
   ylabel('$\theta_x$', 'Interpreter','latex');
   xlabel('Span $x$ [m]', 'Interpreter','latex');

   sgtitle(sprintf('Mode %d  (f = %.2f Hz)', i, frequencies(i)));
end

