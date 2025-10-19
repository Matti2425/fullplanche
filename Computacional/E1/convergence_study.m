%% Convergence Study - CFD Simulation
% Validates numerical solution against analytical solutions for different mesh sizes
% Mesh sizes: N, M = 8, 16, 32, 64 (i.e., 2^3 to 2^6)
%
% This script:
% 1. Runs simulations with increasing mesh refinement
% 2. Computes analytical solutions (3 cases)
% 3. Calculates convergence metrics (L2-norm, L-infinity norm)
% 4. Generates convergence plots
%
clearvars; close all; clc;

%% Setup
mesh_sizes = [8, 16, 32, 64, 2^7, 2^8];  % N and M values: 2^3 to 2^8
num_cases = 3;
num_meshes = length(mesh_sizes);

% Storage for convergence data
convergence_data = struct();
convergence_data.mesh_sizes = mesh_sizes;
convergence_data.L2_error = zeros(num_cases, num_meshes);
convergence_data.Linf_error = zeros(num_cases, num_meshes);
convergence_data.dx = zeros(1, num_meshes);
convergence_data.psi_analytical = cell(num_cases, num_meshes);
convergence_data.psi_numerical = cell(num_cases, num_meshes);
convergence_data.iterations = zeros(num_cases, num_meshes);
convergence_data.convergence_error = zeros(num_cases, num_meshes);

fprintf('=== CONVERGENCE STUDY: CFD SIMULATION ===\n');
fprintf('Mesh sizes (N, M): %s\n', sprintf('%d ', mesh_sizes));
fprintf('Cases: 1=Channel, 2=Static Cylinder, 3=Rotating Cylinder\n\n');

%% Main Convergence Loop
for case_idx = 1:num_cases
    fprintf('---------- CASE %d ----------\n', case_idx);
    
    for mesh_idx = 1:num_meshes
        N = mesh_sizes(mesh_idx);
        M = mesh_sizes(mesh_idx);
        
        fprintf('Mesh: N=%d, M=%d ... ', N, M);
        tic;
        
        %% Read Input Data with Modified Mesh Size
        [Geometry, v_in, p_in, T_in, psi_body_default, gas, Numerical, Init, w, ~] = readInputData();
        Numerical.N = N;
        Numerical.M = M;
        CASE = case_idx;
        
        %% Recompute psi_body for the current case
        if CASE == 1
            psi_body = 0;  % Channel only - no solid
        else
            R_cyl = Geometry.D / 2;
            if CASE == 2
                % Static cylinder: psi_body from uniform flow
                psi_body = v_in * (Geometry.H/2);
            else  % CASE == 3
                % Rotating cylinder: include circulation effect
                psi_body = v_in * (Geometry.H/2) + R_cyl^2 * w * (1/(2*pi));
            end
        end
        
        %% Generate Mesh
        [x_P, y_P] = generateMesh(Geometry, Numerical);
        
        %% Generate Cylinder Geometry and Solid Matrix
        if CASE == 1
            % Channel only - no solid
            S = ones(Numerical.M+2, Numerical.N+2);
            x_cyl = [];
            y_cyl = [];
            theta = [];
        else
            [x_cyl, y_cyl, S, theta] = solidMaterial(Geometry, x_P, y_P, Numerical);
        end
        
        %% Initialize Fields
        [psi, rho] = initParams(Init, Numerical, S);
        
        %% Evaluate Coefficients
        [a_E, a_W, a_N, a_S, b_P, a_P] = evalCoefficients(psi, psi_body, Numerical, Geometry, v_in, x_P, y_P, S, Init);
        
        %% Solve System
        [psi, iter, sol_error] = solveSystem(psi, psi_body, S, a_E, a_W, a_N, a_S, a_P, b_P, v_in, y_P, Numerical, Geometry);
        
        %% Compute Analytical Solution
        R_cyl = Geometry.D / 2;
        psi_analytical = compute_analytical_psi(x_P, y_P, Geometry, v_in, CASE, w, R_cyl);
        
        %% Calculate Error Metrics
        % Only consider fluid domain (outside solid)
        if CASE ~= 1
            error_domain = S == 1;  % Only where fluid exists
        else
            error_domain = ones(size(S));  % All domain for channel case
        end
        
        % Compute L2-norm error
        error_psi = (psi - psi_analytical) .* error_domain;
        L2_error = sqrt(sum(sum(error_psi.^2)) / sum(sum(error_domain)));
        
        % Compute L-infinity norm error
        Linf_error = max(max(abs(error_psi)));
        
        %% Store Results
        convergence_data.L2_error(case_idx, mesh_idx) = L2_error;
        convergence_data.Linf_error(case_idx, mesh_idx) = Linf_error;
        convergence_data.dx(mesh_idx) = Geometry.L / N;
        convergence_data.psi_numerical{case_idx, mesh_idx} = psi;
        convergence_data.psi_analytical{case_idx, mesh_idx} = psi_analytical;
        convergence_data.iterations(case_idx, mesh_idx) = iter;
        convergence_data.convergence_error(case_idx, mesh_idx) = sol_error;
        
        elapsed = toc;
        fprintf('L2-error=%.2e, L∞-error=%.2e, Iter=%d, Time=%.2fs\n', ...
            L2_error, Linf_error, iter, elapsed);
    end
    fprintf('\n');
end

%% Calculate Convergence Rates
fprintf('=== CONVERGENCE RATES ===\n');
fprintf('Case | L2 Rate | Linf Rate | Expected (2nd order)\n');
fprintf('-----|---------|-----------|--------------------\n');

for case_idx = 1:num_cases
    % L2 convergence rate: log(e2/e1) / log(h2/h1)
    if convergence_data.L2_error(case_idx, 1) > 0 && convergence_data.L2_error(case_idx, 2) > 0
        dx1 = convergence_data.dx(1);
        dx2 = convergence_data.dx(2);
        e1 = convergence_data.L2_error(case_idx, 1);
        e2 = convergence_data.L2_error(case_idx, 2);
        
        if e2 > 0
            L2_rate = log(e2/e1) / log(dx2/dx1);
        else
            L2_rate = NaN;
        end
        
        if convergence_data.Linf_error(case_idx, 1) > 0 && convergence_data.Linf_error(case_idx, 2) > 0
            e1_inf = convergence_data.Linf_error(case_idx, 1);
            e2_inf = convergence_data.Linf_error(case_idx, 2);
            Linf_rate = log(e2_inf/e1_inf) / log(dx2/dx1);
        else
            Linf_rate = NaN;
        end
        
        fprintf('%d    |  %.2f   |   %.2f   |         2.0\n', case_idx, L2_rate, Linf_rate);
    end
end

%% Plotting - Convergence Study

% Figure 1: L2-norm convergence
figure('Name', 'Convergence Study - L2-norm', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
hold on;
dx_ref = convergence_data.dx;
for case_idx = 1:num_cases
    loglog(convergence_data.dx, convergence_data.L2_error(case_idx, :), '-', 'LineWidth', 2);
end
% Plot reference lines (1st and 2nd order convergence)
e_ref_1st = convergence_data.L2_error(1, 1) * (dx_ref / dx_ref(1)).^1;
e_ref_2nd = convergence_data.L2_error(1, 1) * (dx_ref / dx_ref(1)).^2;
loglog(dx_ref, e_ref_1st, 'k--', 'LineWidth', 1.5, 'DisplayName', '1st Order');
loglog(dx_ref, e_ref_2nd, 'k:', 'LineWidth', 1.5, 'DisplayName', '2nd Order');
hold off;
xlabel('Grid Spacing (dx)');
ylabel('L2-norm Error');
title('L2-norm Convergence');
legend('Case 1', 'Case 2', 'Case 3', '1st Order', '2nd Order');
grid on;
set(gca, 'XDir', 'reverse');

% Figure 2: L-infinity norm convergence
figure('Name', 'Convergence Study - L∞-norm', 'NumberTitle', 'off', 'Position', [950, 100, 800, 600]);
hold on;
for case_idx = 1:num_cases
    loglog(convergence_data.dx, convergence_data.Linf_error(case_idx, :), '-', 'LineWidth', 2);
end
% Plot reference lines
e_ref_1st_inf = convergence_data.Linf_error(1, 1) * (dx_ref / dx_ref(1)).^1;
e_ref_2nd_inf = convergence_data.Linf_error(1, 1) * (dx_ref / dx_ref(1)).^2;
loglog(dx_ref, e_ref_1st_inf, 'k--', 'LineWidth', 1.5, 'DisplayName', '1st Order');
loglog(dx_ref, e_ref_2nd_inf, 'k:', 'LineWidth', 1.5, 'DisplayName', '2nd Order');
hold off;
xlabel('Grid Spacing (dx)');
ylabel('L∞-norm Error');
title('L∞-norm Convergence');
legend('Case 1', 'Case 2', 'Case 3', '1st Order', '2nd Order');
grid on;
set(gca, 'XDir', 'reverse');

% Figure 3: Iterations vs Mesh Size
figure('Name', 'Convergence Study - Iterations', 'NumberTitle', 'off', 'Position', [100, 750, 800, 600]);
hold on;
for case_idx = 1:num_cases
    semilogy(mesh_sizes, convergence_data.iterations(case_idx, :), '-', 'LineWidth', 2);
end
hold off;
xlabel('Mesh Size (N = M)');
ylabel('Number of Iterations');
title('Solver Iterations vs Mesh Refinement');
legend('Case 1', 'Case 2', 'Case 3');
grid on;
set(gca, 'XScale', 'log');

% Figure 4: Solver Residual vs Mesh Size
figure('Name', 'Convergence Study - Residual', 'NumberTitle', 'off', 'Position', [950, 750, 800, 600]);
hold on;
for case_idx = 1:num_cases
    loglog(mesh_sizes, convergence_data.convergence_error(case_idx, :), '-', 'LineWidth', 2);
end
hold off;
xlabel('Mesh Size (N = M)');
ylabel('Solver Residual');
title('Solver Residual vs Mesh Refinement');
legend('Case 1', 'Case 2', 'Case 3');
grid on;
set(gca, 'XScale', 'log');

%% Print Summary Table
fprintf('\n=== CONVERGENCE SUMMARY TABLE ===\n');
fprintf('Case | Mesh Size | L2-Error   | L∞-Error   | Iterations | Residual\n');
fprintf('-----|-----------|------------|------------|------------|----------\n');
for case_idx = 1:num_cases
    for mesh_idx = 1:num_meshes
        fprintf('%d    | %d x %d     | %.2e | %.2e |    %d      | %.2e\n', ...
            case_idx, mesh_sizes(mesh_idx), mesh_sizes(mesh_idx), ...
            convergence_data.L2_error(case_idx, mesh_idx), ...
            convergence_data.Linf_error(case_idx, mesh_idx), ...
            convergence_data.iterations(case_idx, mesh_idx), ...
            convergence_data.convergence_error(case_idx, mesh_idx));
    end
    fprintf('-----|-----------|------------|------------|------------|----------\n');
end

%% Save Results
save('convergence_study_results.mat', 'convergence_data');
fprintf('\nResults saved to: convergence_study_results.mat\n');
fprintf('Convergence study completed successfully!\n');

%% TOLERANCE SENSITIVITY STUDY
% Study the effect of solver tolerance (Numerical.d) on solution accuracy
fprintf('\n\n=== TOLERANCE SENSITIVITY STUDY ===\n');

tolerance_values = linspace(1e-6, 1e-1, 200)*10e-3;
num_tol = length(tolerance_values);
fixed_mesh_size = 32;  % Use 32x32 mesh for tolerance study

% Storage for tolerance sensitivity data
tolerance_data = struct();
tolerance_data.tolerances = tolerance_values;
tolerance_data.L2_error = zeros(num_cases, num_tol);
tolerance_data.Linf_error = zeros(num_cases, num_tol);
tolerance_data.normalized_L2_error = zeros(num_cases, num_tol);
tolerance_data.normalized_Linf_error = zeros(num_cases, num_tol);
tolerance_data.iterations = zeros(num_cases, num_tol);

fprintf('Fixed mesh size: %d x %d\n', fixed_mesh_size, fixed_mesh_size);
fprintf('Tolerance values: ');
fprintf('%.0e ', tolerance_values);
fprintf('\n\n');

for case_idx = 1:num_cases
    fprintf('---------- CASE %d (Tolerance Study) ----------\n', case_idx);
    
    for tol_idx = 1:num_tol
        tol = tolerance_values(tol_idx);
        
        fprintf('Tolerance: %.0e ... ', tol);
        tic;
        
        %% Read Input Data
        [Geometry, v_in, p_in, T_in, psi_body_default, gas, Numerical, Init, w, ~] = readInputData();
        Numerical.N = fixed_mesh_size;
        Numerical.M = fixed_mesh_size;
        Numerical.d = tol;  % Set tolerance
        CASE = case_idx;
        
        %% Recompute psi_body for the current case
        if CASE == 1
            psi_body = 0;
        else
            R_cyl = Geometry.D / 2;
            if CASE == 2
                psi_body = v_in * (Geometry.H/2);
            else  % CASE == 3
                psi_body = v_in * (Geometry.H/2) + R_cyl^2 * w * (1/(2*pi));
            end
        end
        
        %% Generate Mesh
        [x_P, y_P] = generateMesh(Geometry, Numerical);
        
        %% Generate Cylinder Geometry
        if CASE == 1
            S = ones(Numerical.M+2, Numerical.N+2);
            x_cyl = [];
            y_cyl = [];
            theta = [];
        else
            [x_cyl, y_cyl, S, theta] = solidMaterial(Geometry, x_P, y_P, Numerical);
        end
        
        %% Initialize Fields
        [psi, rho] = initParams(Init, Numerical, S);
        
        %% Evaluate Coefficients
        [a_E, a_W, a_N, a_S, b_P, a_P] = evalCoefficients(psi, psi_body, Numerical, Geometry, v_in, x_P, y_P, S, Init);
        
        %% Solve System
        [psi, iter, sol_error] = solveSystem(psi, psi_body, S, a_E, a_W, a_N, a_S, a_P, b_P, v_in, y_P, Numerical, Geometry);
        
        %% Compute Analytical Solution
        R_cyl = Geometry.D / 2;
        psi_analytical = compute_analytical_psi(x_P, y_P, Geometry, v_in, CASE, w, R_cyl);
        
        %% Calculate Error Metrics
        if CASE ~= 1
            error_domain = S == 1;
        else
            error_domain = ones(size(S));
        end
        
        % Absolute errors
        error_psi = (psi - psi_analytical) .* error_domain;
        L2_error = sqrt(sum(sum(error_psi.^2)) / sum(sum(error_domain)));
        Linf_error = max(max(abs(error_psi)));
        
        % Normalized errors (by maximum analytical psi value)
        psi_max = max(max(abs(psi_analytical(error_domain))));
        if psi_max > 0
            normalized_L2_error = L2_error / psi_max;
            normalized_Linf_error = Linf_error / psi_max;
        else
            normalized_L2_error = L2_error;
            normalized_Linf_error = Linf_error;
        end
        
        %% Store Results
        tolerance_data.L2_error(case_idx, tol_idx) = L2_error;
        tolerance_data.Linf_error(case_idx, tol_idx) = Linf_error;
        tolerance_data.normalized_L2_error(case_idx, tol_idx) = normalized_L2_error;
        tolerance_data.normalized_Linf_error(case_idx, tol_idx) = normalized_Linf_error;
        tolerance_data.iterations(case_idx, tol_idx) = iter;
        
        elapsed = toc;
        fprintf('L2=%.2e, L∞=%.2e, Norm L2=%.2e, Iter=%d, Time=%.2fs\n', ...
            L2_error, Linf_error, normalized_L2_error, iter, elapsed);
    end
    fprintf('\n');
end

%% Plotting - Tolerance Sensitivity Study

% Figure 1: Absolute L2-norm error vs Tolerance
figure('Name', 'Tolerance Study - Absolute L2-norm Error', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
hold on;
for case_idx = 1:num_cases
    semilogx(tolerance_values, tolerance_data.L2_error(case_idx, :), '-', 'LineWidth', 2);
end
hold off;
xlabel('Solver Tolerance (Numerical.d)');
ylabel('L2-norm Error');
title('Absolute L2-norm Error vs Tolerance (Fixed Mesh 32×32)');
legend('Case 1', 'Case 2', 'Case 3');
grid on;

% Figure 2: Absolute L∞-norm error vs Tolerance
figure('Name', 'Tolerance Study - Absolute L∞-norm Error', 'NumberTitle', 'off', 'Position', [950, 100, 800, 600]);
hold on;
for case_idx = 1:num_cases
    semilogx(tolerance_values, tolerance_data.Linf_error(case_idx, :), '-', 'LineWidth', 2);
end
hold off;
xlabel('Solver Tolerance (Numerical.d)');
ylabel('L∞-norm Error');
title('Absolute L∞-norm Error vs Tolerance (Fixed Mesh 32×32)');
legend('Case 1', 'Case 2', 'Case 3');
grid on;

% Figure 3: Normalized L2-norm error vs Tolerance
figure('Name', 'Tolerance Study - Normalized L2-norm Error', 'NumberTitle', 'off', 'Position', [100, 750, 800, 600]);
hold on;
for case_idx = 1:num_cases
    semilogx(tolerance_values, tolerance_data.normalized_L2_error(case_idx, :), '-', 'LineWidth', 2);
end
hold off;
xlabel('Solver Tolerance (Numerical.d)');
ylabel('Normalized L2-norm Error');
title('Normalized L2-norm Error vs Tolerance (Fixed Mesh 32×32)');
legend('Case 1', 'Case 2', 'Case 3');
grid on;

% Figure 4: Normalized L∞-norm error vs Tolerance
figure('Name', 'Tolerance Study - Normalized L∞-norm Error', 'NumberTitle', 'off', 'Position', [950, 750, 800, 600]);
hold on;
for case_idx = 1:num_cases
    semilogx(tolerance_values, tolerance_data.normalized_Linf_error(case_idx, :), '-', 'LineWidth', 2);
end
hold off;
xlabel('Solver Tolerance (Numerical.d)');
ylabel('Normalized L∞-norm Error');
title('Normalized L∞-norm Error vs Tolerance (Fixed Mesh 32×32)');
legend('Case 1', 'Case 2', 'Case 3');
grid on;

%% Print Tolerance Study Summary Table
fprintf('\n=== TOLERANCE SENSITIVITY SUMMARY ===\n');
fprintf('Case | Tolerance  | L2-Error   | L∞-Error   | Norm L2    | Norm L∞    | Iterations\n');
fprintf('-----|------------|------------|------------|------------|------------|------------\n');
for case_idx = 1:num_cases
    for tol_idx = 1:num_tol
        fprintf('%d    | %.0e | %.2e | %.2e | %.2e | %.2e |    %d\n', ...
            case_idx, tolerance_values(tol_idx), ...
            tolerance_data.L2_error(case_idx, tol_idx), ...
            tolerance_data.Linf_error(case_idx, tol_idx), ...
            tolerance_data.normalized_L2_error(case_idx, tol_idx), ...
            tolerance_data.normalized_Linf_error(case_idx, tol_idx), ...
            tolerance_data.iterations(case_idx, tol_idx));
    end
    fprintf('-----|------------|------------|------------|------------|------------|------------\n');
end

%% Save Tolerance Study Results
save('tolerance_sensitivity_results.mat', 'tolerance_data');
fprintf('\nTolerance study results saved to: tolerance_sensitivity_results.mat\n');