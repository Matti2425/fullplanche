%% PART II - SHELL MODELLING
% Gerard Martínez Vilella
% David Rodríguez Gonzalez
% Useful commands to initialize script
clear
close all
clc

%% DATA
problem_type = 'bending'; % bending or torsion
Nm = 6; % Number of modes to extract

% Dimensions and Material Properties
E = 75e9; % Pa
v_poisson = 0.3;
rho = 2200; % kg/m3

% Load definition
F_load = 1; % N (Force)
T_load = 1; % Nm (Torque)

% Geometric definitions
c = 0.85; b = 3;
y0 = 0.924; %m
y1 = 0.125; %m
y2 = 0.475; %m
yc = 0.285; %m % Shear center approx
h1 = 35e-3; %m
h2 = 30e-3; %m
h3 = 25e-3; %m

%% PREPROCESS
% Load mesh data
% Ensure 'shell.mat' exists in your path
try
    load('shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');
catch
    error('File shell.mat not found. Please ensure the mesh file is in the directory.');
end

N_elements = size(Tn,1);
N_nodes = size(xn,1);
Ndof = 6 * N_nodes;

% --- 1. Define Boundary Conditions (Up) ---
% Fix all 6 DOFs at the root nodes
% Vectorized approach for efficiency
Up = zeros(length(indRoot)*6, 3);
idx = 1;
for i = 1:length(indRoot)
    for dof = 1:6
        Up(idx,:) = [0, indRoot(i), dof]; % Displacement, Node, DOF
        idx = idx + 1;
    end
end

% --- 2. Define Loads (Fe) ---
% Calculate load distribution based on geometry
Fa_shear = F_load * ((y2-yc)/(y2-y1));
Fb_shear = F_load * ((yc-y1)/(y2-y1));
Fa_torque = -T_load / (y2-y1);
Fb_torque = T_load / (y2-y1);

switch(problem_type)
    case 'bending' % Bending
        % Force in Z (DOF 3)
        Fe = [Fa_shear, indPointA, 3;
              Fb_shear, indPointB, 3];
    case 'torsion' % Torsion
        % Force in Z (DOF 3) creating a couple
        Fe = [Fa_torque, indPointA, 3;
              Fb_torque, indPointB, 3];
end

% Placeholder for distributed/body forces (empty as per reference)
Pe = zeros(0,3); 
Be = zeros(0,3);

% Property Vectors based on material index Tm
% Assuming Tm maps to: 1=h1, 2=h2, 3=h3
h_vals = [h1; h2; h3];
h_vector = h_vals(Tm);
E_vector = repmat(E, N_elements, 1);
v_vector = repmat(v_poisson, N_elements, 1);
rho_vector = repmat(rho, N_elements, 1);

%% SOLVER
if 1 % Set to 1 to (re)compute, 0 to load
    fprintf('--- Starting Matrix Assembly ---\n');
    
    K = sparse(Ndof, Ndof);
    M = sparse(Ndof, Ndof);
    
    % Pre-computations for Gauss integration (MITC4 / Q4 Shell logic)
    % Standard linear shape function derivatives for Q4
    a_coef = [-1 1 1 -1];
    b_coef = [-1 -1 1 1];
    
    % Integration points for Mass and Bending (2x2 Gauss)
    xi4 = [-1 1 1 -1]/sqrt(3);
    nu4 = [-1 -1 1 1]/sqrt(3);
    w4 = [1 1 1 1];
    
    for e = 1:N_elements
        % --- 1. Local Coordinate System (R) ---
        % Nodes: 1-2-3-4
        x_e = xn(Tn(e,:), :); % 4x3 coordinates
        
        % Vectors for local axes
        v13 = x_e(3,:) - x_e(1,:);
        v24 = x_e(4,:) - x_e(2,:);
        Normal = cross(v13, v24); 
        k_dir = (Normal / norm(Normal))'; % Normal to element (local z)
        
        % Local x (i_dir) approx aligned with 1-2 or midpoints
        v_bisect = (x_e(2,:) + x_e(3,:) - x_e(4,:) - x_e(1,:)) / 2;
        i_dir = (v_bisect / norm(v_bisect))';
        j_dir = cross(k_dir, i_dir); % Local y
        
        % Rotation Matrix (3x3) -> Expanded to (24x24)
        Rp = [i_dir, j_dir, k_dir]; % 3x3 global to local
        % Full element rotation matrix (blkdiag logic manually implemented for speed)
        % Note: The reference code assumes 5 local DOFs + drilling, mapped to 6 global
        % But actually constructs a 24x24 matrix (4 nodes * 6 dofs)
        
        % Usually R maps Global U to Local U. u_local = R' * u_global ?
        % The reference defines Rp such that [u_local] = Rp' * [u_global]
        % Let's build the 3x3 block:
        R_block = blkdiag(Rp', Rp'); % 6x6 block (Rotation for u and theta)
        % Actually, usually we rotate the stiffness: K_global = R * K_local * R'
        % The reference creates a matrix R of size 20x24. 
        % 5 local DOFs (u,v,w, th_x, th_y) per node -> 20 rows
        % 6 global DOFs per node -> 24 columns.
        % This implies the local formulation ignores the "drilling" rotation initially.
        
        % Reconstructing Reference Logic exactly:
        R_node = zeros(5,6);
        R_node(1:3, 1:3) = Rp'; % Rotates displacements
        R_node(4:5, 4:6) = Rp(:,1:2)'; % Rotates rotations (takes only local x,y components)
        
        R_ele = zeros(20,24);
        for n = 1:4
            R_ele( (n-1)*5+(1:5), (n-1)*6+(1:6) ) = R_node;
        end
        
        % --- 2. Stiffness Matrices ---
        
        % A) Shear Stiffness (Ks) - Reduced integration / MITC assumptions
        % Center point (xi=0, nu=0) evaluation for shear (avoid locking)
        J_center = zeros(2,2);
        for n=1:4
            dNi_dxi = [a_coef(n); b_coef(n)]/4;
            J_center = J_center + dNi_dxi * x_e(n,:) * [i_dir, j_dir];
        end
        detJ = det(J_center);
        invJ = inv(J_center);
        
        Bs = zeros(2, 20); % 2 strains (gamma_xz, gamma_yz), 20 DOFs
        for n=1:4
            % Derivatives of N mapped to physical x,y
            dN_phys = invJ * [a_coef(n); b_coef(n)]/4; 
            N_val = 1/4; % at center
            
            % Shear B-matrix part for node n
            % Local DOFs: [u, v, w, th_x, th_y]
            % Strain-disp: gamma_xz = dw/dx + th_y ... gamma_yz = dw/dy - th_x
            % (Check specific element formulation signs)
            % Reference:
            % Row 1 (xz): N,x * w + N * th_y 
            % Row 2 (yz): N,y * w - N * th_x
            
            col_idx = (n-1)*5 + (1:5);
            %              u v  w           th_x    th_y
            Bs(1, col_idx) = [0 0 dN_phys(1)  0       N_val];
            Bs(2, col_idx) = [0 0 dN_phys(2) -N_val   0    ];
        end
        
        Cs = eye(2) * 5/6 * E_vector(e)*h_vector(e) / (2*(1+v_vector(e))); % Shear modulus * Area * k
        Ks = (4 * detJ) * (R_ele' * Bs' * Cs * Bs * R_ele); % Area = 4*detJ roughly
        
        % B) Membrane (Km) and Bending (Kb) + Mass (Me) - Full integration
        Km = sparse(24,24);
        Kb = sparse(24,24);
        Me = sparse(24,24);
        
        C_membrane = E_vector(e)*h_vector(e)/(1-v_vector(e)^2) * ...
                     [1 v_vector(e) 0; v_vector(e) 1 0; 0 0 (1-v_vector(e))/2];
                 
        C_bending  = E_vector(e)*h_vector(e)^3/(12*(1-v_vector(e)^2)) * ...
                     [1 v_vector(e) 0; v_vector(e) 1 0; 0 0 (1-v_vector(e))/2];
                 
        rho_mat_base = rho_vector(e)*h_vector(e)*diag([1 1 1 h_vector(e)^2/12 h_vector(e)^2/12]);

        for k = 1:4 % Loop over 4 Gauss points
            % Shape functions and derivatives at Gauss point k
            xi = xi4(k); nu = nu4(k);
            
            J = zeros(2,2);
            N = zeros(1,4);
            dN_dxi = zeros(2,4);
            
            for n=1:4
                val = (1 + a_coef(n)*xi)*(1 + b_coef(n)*nu)/4;
                dval_dxi = a_coef(n)*(1 + b_coef(n)*nu)/4;
                dval_dnu = b_coef(n)*(1 + a_coef(n)*xi)/4;
                
                N(n) = val;
                dN_dxi(:,n) = [dval_dxi; dval_dnu];
                
                J = J + dN_dxi(:,n) * x_e(n,:) * [i_dir, j_dir];
            end
            
            detJ_k = det(J);
            dN_phys = J \ dN_dxi; % 2x4
            
            Bm = zeros(3, 20); % Membrane strains
            Bb = zeros(3, 20); % Curvatures
            N_mat = zeros(5, 20); % Mass interpolation
            
            for n=1:4
                col = (n-1)*5 + (1:5);
                % Membrane [du/dx; dv/dy; du/dy+dv/dx]
                Bm(1, col) = [dN_phys(1,n) 0 0 0 0];
                Bm(2, col) = [0 dN_phys(2,n) 0 0 0];
                Bm(3, col) = [dN_phys(2,n) dN_phys(1,n) 0 0 0];
                
                % Bending [dth_y/dx; -dth_x/dy; dth_y/dy - dth_x/dx]
                % Note: Assuming th_x is rot about x, th_y about y.
                % Reference logic check:
                % Row 1: d(th_y)/dx -> col 5
                % Row 2: -d(th_x)/dy -> col 4
                Bb(1, col) = [0 0 0  0             dN_phys(1,n)];
                Bb(2, col) = [0 0 0 -dN_phys(2,n)  0           ];
                Bb(3, col) = [0 0 0 -dN_phys(1,n)  dN_phys(2,n)];
                
                % Mass matrix shape function
                N_mat(:, col) = N(n) * eye(5);
            end
            
            dV = w4(k) * detJ_k;
            
            Km = Km + dV * (R_ele' * Bm' * C_membrane * Bm * R_ele);
            Kb = Kb + dV * (R_ele' * Bb' * C_bending * Bb * R_ele);
            Me = Me + dV * (R_ele' * N_mat' * rho_mat_base * N_mat * R_ele);
        end
        
        % Assemble Element
        K_el = Ks + Km + Kb;
        
        % Global Indices
        Idof = [];
        for n=1:4
            Idof = [Idof, 6*(Tn(e,n)-1) + (1:6)];
        end
        
        K(Idof, Idof) = K(Idof, Idof) + K_el;
        M(Idof, Idof) = M(Idof, Idof) + Me;
    end
    
    % --- 3. Drilling Rotation Stabilization (Kr) ---
    % Add artificial stiffness for the 6th DOF (Rotation about normal Z)
    % otherwise the matrix is singular for flat shells.
    fprintf('Applying Drilling Rotation Stiffness...\n');
    % Calculate node normals
    Normal_Node = zeros(N_nodes, 3);
    for e=1:N_elements
        v13 = xn(Tn(e,3),:) - xn(Tn(e,1),:);
        v24 = xn(Tn(e,4),:) - xn(Tn(e,2),:);
        n_ele = cross(v13, v24);
        n_ele = n_ele / norm(n_ele);
        for i=1:4
            Normal_Node(Tn(e,i), :) = Normal_Node(Tn(e,i), :) + n_ele;
        end
    end
    
    Kr = sparse(Ndof, Ndof);
    alpha_coeff = 1e-3 * E * h1; % Small stiffness coefficient
    
    for e=1:N_elements
        % Re-calculate element normal
        v13 = xn(Tn(e,3),:) - xn(Tn(e,1),:);
        v24 = xn(Tn(e,4),:) - xn(Tn(e,2),:);
        n_ele = cross(v13, v24);
        Area = norm(n_ele)/2;
        n_ele = n_ele / norm(n_ele);
        
        for i=1:4
            node = Tn(e,i);
            n_node = Normal_Node(node, :) / norm(Normal_Node(node, :));
            
            % Check alignment. If element normal is close to node normal (smooth area)
            % add stiffness to drilling dof.
            alignment = dot(n_ele, n_node);
            
            % If reasonably flat (alignment ~ 1), add penalty
            if abs(alignment) > 0.9
                % Project stiffness onto global DOFs [th_x, th_y, th_z]
                % Local drilling is rotation about n_ele
                % Stiffness matrix part: k * (n_ele * n_ele')
                K_rot_penalty = (alpha_coeff * Area) * (n_ele' * n_ele);
                
                idx_rot = 6*(node-1) + (4:6);
                Kr(idx_rot, idx_rot) = Kr(idx_rot, idx_rot) + K_rot_penalty;
            end
        end
    end
    K = K + Kr;

    save('shell_matrices.mat','K','M'); 
    fprintf('System Matrices Computed and Saved.\n');
else
    load('shell_matrices.mat','K','M');
    fprintf('System Matrices Loaded.\n');
end

% --- 4. Apply Boundary Conditions and Solve ---
fprintf('Solving Static System...\n');

% Construct Global Force Vector
F_global = zeros(Ndof, 1);
for i = 1:size(Fe,1)
    node = Fe(i,2);
    dof = Fe(i,3);
    val = Fe(i,1);
    idx = 6*(node-1) + dof;
    F_global(idx) = F_global(idx) + val;
end

% Partition
prescribed_dofs = 6*(Up(:,2)-1) + Up(:,3);
free_dofs = setdiff(1:Ndof, prescribed_dofs);

% Solve U
U = zeros(Ndof, 1);
U(prescribed_dofs) = Up(:,1); % Non-zero prescribed displacements if any

% F_free_eff = F_free - K_fp * U_p
F_eff = F_global(free_dofs) - K(free_dofs, prescribed_dofs) * U(prescribed_dofs);
U(free_dofs) = K(free_dofs, free_dofs) \ F_eff;

fprintf('Static Solution Complete.\n');

% --- 5. Modal Analysis ---
fprintf('Performing Modal Analysis...\n');
n_modes_calc = min(Nm + 4, length(free_dofs)); % Calculate a few more for stability
[V_red, D_red] = eigs(K(free_dofs, free_dofs), M(free_dofs, free_dofs), n_modes_calc, 'sm');

[frequencies, sort_idx] = sort(sqrt(diag(D_red))/(2*pi));
V_red = V_red(:, sort_idx);

% Reconstruct full eigenvectors
Phi = zeros(Ndof, n_modes_calc);
Phi(free_dofs, :) = V_red;
freq = frequencies(1:Nm); % Keep requested number

fprintf('Natural Frequencies (Hz):\n');
disp(freq);

%% POSTPROCESS
save('shell_results.mat');

% --- 1. Spanwise Calculation (Specific Request) ---
% Extract results along spars
n_pts = size(indSpar1,1);
res_z = zeros(n_pts,1);
res_theta = zeros(n_pts,1);
res_x = xn(indSpar1, 1);

for i=1:n_pts
    node1 = indSpar1(i);
    node2 = indSpar2(i);
    
    idx1_z = 6*(node1-1) + 3;
    idx2_z = 6*(node2-1) + 3;
    
    uz1 = U(idx1_z);
    uz2 = U(idx2_z);
    
    % Twist (angle between spars)
    % theta ~ (uz_rear - uz_front) / distance_y
    % Or strictly from rotation DOFs if reliable. Let's use displacements as per reference logic.
    dist_y = y2 - y1;
    theta_val = (uz2 - uz1) / dist_y;
    
    % Vertical deflection at Shear Center (extrapolated)
    % z_sc = z1 + theta * (y_sc - y1)
    z_sc = uz1 + theta_val * (yc - y1);
    
    res_z(i) = z_sc;
    res_theta(i) = theta_val;
end

% Plot Spanwise Results
figure('Name','Spanwise Results');
subplot(2,1,1);
plot(res_x, res_z, 'b-o', 'LineWidth', 1.5);
title('Vertical Deflection (u_z) at Shear Center');
xlabel('Span (m)'); ylabel('u_z (m)'); grid on;
legend(sprintf('Load Case %d', problem_type));

subplot(2,1,2);
plot(res_x, res_theta, 'r-s', 'LineWidth', 1.5);
title('Twist Angle (\theta_x)');
xlabel('Span (m)'); ylabel('\theta_x (rad)'); grid on;

% --- 2. 3D Visualizations using requested functions ---

% Plot Deformed Shape
scale_factor = 0.1 * max(abs(xn(:,1))) / max(abs(U(3:6:end))); % Auto-scale approx
if isnan(scale_factor), scale_factor=1; end
fprintf('Plotting deformed shape with scale factor: %.2f\n', scale_factor);
plotDeformed('Shell Static Deformation', xn, Tn, U, scale_factor);

% Plot Modes
imodes = 1:min(6, Nm);
plotModes('Modes Shapes', Phi, freq, imodes);


% MODE ANALYSIS
%% --- Modal postprocess (Part II requirement): ubar_y, ubar_z, thetabar_x vs span ---

n_pts = length(indSpar1);
x_span = xn(indSpar1,1);
dist_y = (y2 - y1);

for k = 1:Nm   % or for k = 1:n_modes_calc if you prefer
    % Modal vector (full size Ndof)
    phi = Phi(:,k);

    ubar_y = zeros(n_pts,1);
    ubar_z = zeros(n_pts,1);
    thetabar_x = zeros(n_pts,1);

    for i = 1:n_pts
        n1 = indSpar1(i);   % front spar node index at station i
        n2 = indSpar2(i);   % rear spar node index at station i

        uy1 = phi(6*(n1-1) + 2);
        uy2 = phi(6*(n2-1) + 2);

        uz1 = phi(6*(n1-1) + 3);
        uz2 = phi(6*(n2-1) + 3);

        % Averages / twist (statement formulas)
        ubar_y(i) = 0.5*(uy1 + uy2);
        thetabar_x(i) = (uz2 - uz1) / dist_y;
        ubar_z(i) = uz1 + thetabar_x(i)*(yc - y1);
    end

    % Normalize for visualization (optional, same reason as Part I)
    scale = max(abs([ubar_y; ubar_z; thetabar_x])) + eps;
    ubar_y = ubar_y/scale;
    ubar_z = ubar_z/scale;
    thetabar_x = thetabar_x/scale;

    figure;
    subplot(3,1,1); plot(x_span, ubar_y, 'LineWidth', 1.5); grid on;
    ylabel('$\bar{u}_y$', 'Interpreter','latex', 'FontSize', 16);

    subplot(3,1,2); plot(x_span, ubar_z, 'LineWidth', 1.5); grid on;
    ylabel('$\bar{u}_z$', 'Interpreter','latex', 'FontSize', 16);

    subplot(3,1,3); plot(x_span, thetabar_x, 'LineWidth', 1.5); grid on;
    ylabel('$\bar{\theta}_x$', 'Interpreter','latex', 'FontSize', 16);
    xlabel('Span $x$ [m]', 'Interpreter','latex', 'FontSize', 16);

    sgtitle(sprintf('Mode %d (f = %.2f Hz)', k, freq(k)));
end


% MODE ANALYSYS


%% MODIFIED PLOT MODES AND PLOTDEFORMED

function plotDeformed(plotName, xn, Tn, u, scale)
    % xn : [Nnodes x 3]
    % Tn : [Nelem x 4]
    % u  : [Ndof x 1]
    
    figure('Name', plotName, 'NumberTitle', 'off', 'Color', 'w');
    
    % Extract displacements
    ux = u(1:6:end);
    uy = u(2:6:end);
    uz = u(3:6:end);
    
    % Deformed coordinates
    xn_def = xn + scale * [ux, uy, uz];
    
    % Plot undeformed (wireframe, transparent)
    patch('Vertices', xn, 'Faces', Tn, ...
          'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7], 'LineStyle', ':');
    hold on;
    
    % Plot deformed (colored by Z displacement)
    % CData is used for color mapping. We map Z-displacement to color.
    % We need one color value per vertex or per face.
    cdata = uz; 
    
    p = patch('Vertices', xn_def, 'Faces', Tn, ...
          'FaceVertexCData', cdata, ...
          'FaceColor', 'interp', 'EdgeColor', 'k', 'FaceAlpha', 0.9);
      
    title(sprintf('%s (Scale: %.1f)', plotName, scale));
    axis equal; axis vis3d;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    colorbar; colormap jet;
    view(3); grid on;
end

function plotModes(plotName, Phi, freq, imodes)
    % Phi    : [Ndof x Nmodes]
    % freq   : [Nmodes x 1]
    % imodes : Vector of indices to plot
    
    % Load mesh data locally or pass it in. 
    % Since functions define scope, we assume xn, Tn are available via evalin 
    % or we need to pass them. For this template, I will assume variables exist 
    % in workspace (MATLAB script behavior) or use evalin. 
    % Better practice: Pass xn/Tn as arguments. I will fetch from base workspace.
    
    xn = evalin('base', 'xn');
    Tn = evalin('base', 'Tn');
    
    figure('Name', plotName, 'NumberTitle', 'off', 'Color', 'w');
    
    num_plots = length(imodes);
    rows = ceil(sqrt(num_plots));
    cols = ceil(num_plots/rows);
    
    for k = 1:num_plots
        mode_idx = imodes(k);
        subplot(rows, cols, k);
        
        mode_shape = Phi(:, mode_idx);
        
        % Normalize mode for plotting
        max_disp = max(abs(mode_shape(1:6:end).^2 + mode_shape(2:6:end).^2 + mode_shape(3:6:end).^2));
        scale = 0.1 * max(max(xn)) / sqrt(max_disp);
        
        % Extract deformation
        dx = mode_shape(1:6:end);
        dy = mode_shape(2:6:end);
        dz = mode_shape(3:6:end);
        
        xn_def = xn + scale * [dx, dy, dz];
        
        % Plot Undeformed Shape (Wireframe)
        patch('Vertices', xn, 'Faces', Tn, 'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7], 'LineStyle', ':');
        hold on;
        
        % Plot Deformed Shape (Colored by Z displacement)
        patch('Vertices', xn_def, 'Faces', Tn, 'FaceVertexCData', dz, ...
              'FaceColor', 'interp', 'EdgeColor', 'k', 'FaceAlpha', 0.8);
          
        % Set the title (static and positioned outside the plot)
        title(sprintf('Mode %d: %.2f Hz', mode_idx, freq(mode_idx)), ...
              'Position', [1, -0.6, -0.5], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        axis equal; 
        axis vis3d; 
        axis off;
        view(3);
    end
    sgtitle(plotName);
    colormap jet;
end