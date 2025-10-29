%% CFD Simulation - Generic Convection-Diffusion Transport Equation
% E2 Implementation
% Assignment 2: Solving the generic convection-diffusion equation with various schemes
clearvars; close all;

%% CONFIGURATION
% Problem type:
% diag = false: Smith-Hutton test case
% diag = true: Diagonal flow test case
diag = false;

% Override for diagnostic mode
if diag
    case_study = "UDS";  % Use Upwind for diagonal case
else
    % Case study selection:
    % UDS (Upwind Differencing Scheme)
    % CDS (Central Differencing Scheme)
    % HDS (Hybrid Differencing Scheme)
    % EDS (Exponential Differencing Scheme)
    % PLDS (Power Law Differencing Scheme)
    case_study = "PLDS";
end


%% READ INPUT DATA
fprintf('=== CFD Convection-Diffusion Solver ===\n');
fprintf('Case study: %s\n', case_study);
fprintf('Problem type: %s\n\n', iif(diag, 'Diagonal flow', 'Smith-Hutton'));

[Geometry, Fluid, Numerical, BoundaryConditions] = setupSimulation(diag);

%% GENERATE MESH
fprintf('Generating mesh...\n');
[nod_x, nod_y, center_x, center_y, sx_w, sx_e, sx_p, sy_s, sy_n, sy_p] = generateMesh(Geometry, Numerical);

% Cell dimensions
A_x = Geometry.L / Numerical.N;
A_y = Geometry.H / Numerical.M;

% Create meshgrid for visualization and velocity computation
x = nod_x - Geometry.H;
[surf_x, surf_y] = meshgrid(x, fliplr(nod_y));
surf_x = surf_x';
surf_y = fliplr(surf_y');

fprintf('Mesh: %d x %d cells\n', Numerical.N, Numerical.M);
fprintf('Cell size: dx = %.4f m, dy = %.4f m\n\n', A_x, A_y);

%% COMPUTE VELOCITY FIELD
fprintf('Computing velocity field...\n');
[v_x, v_y] = computeVelocityField(surf_x, surf_y, Fluid, diag);

%% SOLVER LOOP FOR DIFFERENT GAMMA VALUES
phi_def = zeros(Numerical.N+2, Numerical.M+2, length(Fluid.gamma));

for l = 1:length(Fluid.gamma)
    fprintf('--- Solving for rho/Gamma = %.0e ---\n', Fluid.GAMMA(l));
    
    % Initialize phi field
    phi_in = ones(Numerical.N+2, Numerical.M+2) * BoundaryConditions.phi_inlet;
    phi_prev = phi_in;
    
    % Evaluate coefficients
    fprintf('Evaluating coefficients...\n');
    [a_n, a_s, a_e, a_w, a_p, a_p_base, b_p] = computeCoefficients(v_x, v_y, Numerical, Fluid, ...
                                                                  sx_w, sx_e, sx_p, sy_s, sy_n, sy_p, ...
                                                                  Fluid.gamma(l), case_study, A_x, A_y);
    
    % Apply boundary conditions
    fprintf('Applying boundary conditions...\n');
    [a_n, a_s, a_e, a_w, a_p, b_p] = applyBoundaryConditions(a_n, a_s, a_e, a_w, a_p, a_p_base, b_p, ...
                                                              surf_x, Numerical, diag, BoundaryConditions);
    
    % Solve system
    fprintf('Solving system...\n');
    [phi, iter_outer, iter_inner] = solveSystem(a_n, a_s, a_e, a_w, a_p, a_p_base, b_p, ...
                                                 phi_in, phi_prev, Numerical);
    
    fprintf('Converged: Outer iterations = %d, Total inner iterations = %d\n\n', iter_outer, iter_inner);
    
    % Store solution
    phi_def(:, :, l) = phi;
end

%% PLOT RESULTS
fprintf('Generating plots...\n');
plotResults(surf_x, surf_y, phi_def, nod_x, nod_y, Geometry, Fluid, diag, case_study);

fprintf('\n=== Simulation completed successfully! ===\n');

%% HELPER FUNCTION
function out = iif(condition, true_val, false_val)
    % Inline if-else function
    if condition
        out = true_val;
    else
        out = false_val;
    end
end
