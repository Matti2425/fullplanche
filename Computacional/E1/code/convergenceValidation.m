%% Numerical Convergence Analysis and Validation Study
% Validates numerical solutions against analytical solutions with mesh refinement
% Tests convergence rates and tolerance effects
%
% Tests convergence for case 1 (channel) across different mesh resolutions
%
clearvars; close all; clc;

%% Configuration Phase
gridSizes = [8, 16, 32, 64, 128];  % Mesh refinement levels: 2^3 to 2^8
numCases = 3;
numGrids = length(gridSizes);

% Data storage structure
convergenceMetrics = struct();
convergenceMetrics.grid_sizes = gridSizes;
convergenceMetrics.error_L2 = zeros(numCases, numGrids);
convergenceMetrics.error_Linf = zeros(numCases, numGrids);
convergenceMetrics.spacing_dx = zeros(1, numGrids);
convergenceMetrics.psi_analytical_ref = cell(numCases, numGrids);
convergenceMetrics.psi_numerical_ref = cell(numCases, numGrids);
convergenceMetrics.solver_iterations = zeros(numCases, numGrids);
convergenceMetrics.solver_residual = zeros(numCases, numGrids);

fprintf('========= CONVERGENCE ANALYSIS FRAMEWORK =========\n');
fprintf('Grid sizes (N, M): %s\n', sprintf('%d ', gridSizes));
fprintf('Cases: 1=Channel, 2=Static Cylinder, 3=Rotating Cylinder\n\n');

%% Main Convergence Refinement Loop
for caseIdx = 1:numCases
    fprintf('---------- CASE %d CONVERGENCE ----------\n', caseIdx);

    for gridIdx = 1:numGrids
        N = gridSizes(gridIdx);
        M = gridSizes(gridIdx);

        fprintf('  Mesh: %d×%d ... ', N, M);
        timeStart = tic;

        % Initialize Simulation Parameters
        [domainGeo, flowInlet, flowBoundary, fluidMaterial, numericalCfg, ...
            initialState, cylinderSpeed] = setupSimulation(caseIdx);

        numericalCfg.nx = N;
        numericalCfg.ny = M;
        flowCase = caseIdx;

        % Recompute boundary condition for current case
        if flowCase == 1
            psiBoundary = 0;  % Open channel
        else
            r_cyl = domainGeo.diameter / 2;
            if flowCase == 2
                psiBoundary = flowInlet.velocity * (domainGeo.height/2);
            else  % flowCase == 3
                psiBoundary = flowInlet.velocity * (domainGeo.height/2) + ...
                             r_cyl^2 * cylinderSpeed * (1/(2*pi));
            end
        end

        % Mesh Generation
        [meshX, meshY] = createComputationalGrid(domainGeo, numericalCfg);

        % Obstacle Definition
        if flowCase == 1
            % Channel without obstacle
            regionFluid = ones(numericalCfg.ny+2, numericalCfg.nx+2);
            obstacleX = [];
            obstacleY = [];
            anglePos = [];
        else
            [obstacleX, obstacleY, regionFluid, anglePos] = ...
                defineObstacles(domainGeo, meshX, meshY, numericalCfg);
        end

        % Field Initialization
        [psiField, rhoField] = initializeFields(initialState, numericalCfg, regionFluid);

        % Discretization Matrix Construction
        [coeff_e, coeff_w, coeff_n, coeff_s, coeff_c, rhs_vec] = ...
            constructMatrix(psiField, psiBoundary, numericalCfg, domainGeo, ...
                           flowInlet.velocity, meshX, meshY, regionFluid, initialState);

        % Solve Stream Function
        [psiField, numIterations, finalResidual] = iterativeSolver(psiField, psiBoundary, ...
                                                                   regionFluid, coeff_e, coeff_w, ...
                                                                   coeff_n, coeff_s, coeff_c, rhs_vec, ...
                                                                   flowInlet.velocity, meshY, numericalCfg, domainGeo);

        % Reference Analytical Solution
        radiusObs = domainGeo.diameter / 2;
        psiAnalyticalRef = streamFunctionAnalytical(meshX, meshY, domainGeo, ...
                                                    flowInlet.velocity, flowCase, ...
                                                    cylinderSpeed, radiusObs);

        % Error Metric Calculation
        if flowCase ~= 1
            validRegion = regionFluid == 1;  % Only compute error in fluid domain
        else
            validRegion = ones(size(regionFluid));  % All domain for channel
        end

        % L2-norm Error
        errorDomain = (psiField - psiAnalyticalRef) .* validRegion;
        errorL2 = sqrt(sum(sum(errorDomain.^2)) / sum(sum(validRegion)));

        % L-infinity Norm Error
        errorLinf = max(max(abs(errorDomain)));

        % Store Results
        convergenceMetrics.error_L2(caseIdx, gridIdx) = errorL2;
        convergenceMetrics.error_Linf(caseIdx, gridIdx) = errorLinf;
        convergenceMetrics.spacing_dx(gridIdx) = domainGeo.length / N;
        convergenceMetrics.psi_numerical_ref{caseIdx, gridIdx} = psiField;
        convergenceMetrics.psi_analytical_ref{caseIdx, gridIdx} = psiAnalyticalRef;
        convergenceMetrics.solver_iterations(caseIdx, gridIdx) = numIterations;
        convergenceMetrics.solver_residual(caseIdx, gridIdx) = finalResidual;

        timeElapsed = toc(timeStart);
        fprintf('L2-error=%.2e, L∞-error=%.2e, Iter=%d, Time=%.2fs\n', ...
            errorL2, errorLinf, numIterations, timeElapsed);
    end
    fprintf('\n');
end

%% Convergence Rate Computation
fprintf('=== CONVERGENCE RATE ANALYSIS ===\n');
fprintf('Case | L2 Order | L∞ Order | Expected (2nd-order scheme)\n');
fprintf('-----|----------|----------|---------------------------\n');

for caseIdx = 1:numCases
    if convergenceMetrics.error_L2(caseIdx, 1) > 0 && convergenceMetrics.error_L2(caseIdx, 2) > 0
        dx1 = convergenceMetrics.spacing_dx(1);
        dx2 = convergenceMetrics.spacing_dx(2);
        e1 = convergenceMetrics.error_L2(caseIdx, 1);
        e2 = convergenceMetrics.error_L2(caseIdx, 2);

        if e2 > 0
            orderL2 = log(e2/e1) / log(dx2/dx1);
        else
            orderL2 = NaN;
        end

        if convergenceMetrics.error_Linf(caseIdx, 1) > 0 && convergenceMetrics.error_Linf(caseIdx, 2) > 0
            e1_inf = convergenceMetrics.error_Linf(caseIdx, 1);
            e2_inf = convergenceMetrics.error_Linf(caseIdx, 2);
            orderLinf = log(e2_inf/e1_inf) / log(dx2/dx1);
        else
            orderLinf = NaN;
        end

        fprintf('%d    |  %.2f    |   %.2f   |           2.0\n', caseIdx, orderL2, orderLinf);
    end
end

%% Visualization: Convergence Performance
figPos = [100, 100, 900, 700];

% Error vs Grid Spacing (L2-norm)
figure('Name', 'L2-norm Convergence (Case 1)', 'NumberTitle', 'off', 'Position', figPos);
hold on;
dxRef = convergenceMetrics.spacing_dx;
loglog(convergenceMetrics.spacing_dx, convergenceMetrics.error_L2(1, :), ...
       'o-', 'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', 'Numerical');

% Reference convergence lines
eRef1st = convergenceMetrics.error_L2(1, 1) * (dxRef / dxRef(1)).^1;
eRef2nd = convergenceMetrics.error_L2(1, 1) * (dxRef / dxRef(1)).^2;
loglog(dxRef, eRef1st, 'k--', 'LineWidth', 1.5, 'DisplayName', '1st Order');
loglog(dxRef, eRef2nd, 'k:', 'LineWidth', 1.5, 'DisplayName', '2nd Order');
hold off;
xlabel('Grid Spacing (dx)');
ylabel('L2-norm Error');
title('L2-norm Convergence vs Mesh Refinement (Case 1 - Channel Flow)');
legend('FontSize', 11, 'Location', 'best');
grid on;
set(gca, 'XDir', 'reverse');

% Error vs Grid Spacing (L-infinity norm)
figure('Name', 'L∞-norm Convergence (Case 1)', 'NumberTitle', 'off', 'Position', figPos);
hold on;
loglog(convergenceMetrics.spacing_dx, convergenceMetrics.error_Linf(1, :), ...
       's-', 'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', 'Numerical');

eRef1stInf = convergenceMetrics.error_Linf(1, 1) * (dxRef / dxRef(1)).^1;
eRef2ndInf = convergenceMetrics.error_Linf(1, 1) * (dxRef / dxRef(1)).^2;
loglog(dxRef, eRef1stInf, 'k--', 'LineWidth', 1.5, 'DisplayName', '1st Order');
loglog(dxRef, eRef2ndInf, 'k:', 'LineWidth', 1.5, 'DisplayName', '2nd Order');
hold off;
xlabel('Grid Spacing (dx)');
ylabel('L∞-norm Error');
title('L∞-norm Convergence vs Mesh Refinement (Case 1 - Channel Flow)');
legend('FontSize', 11, 'Location', 'best');
grid on;
set(gca, 'XDir', 'reverse');

% Solver Iterations vs Mesh Size
figure('Name', 'Iterations vs Mesh Refinement (Case 1)', 'NumberTitle', 'off', 'Position', figPos);
semilogy(gridSizes, convergenceMetrics.solver_iterations(1, :), ...
         'o-', 'LineWidth', 2.5, 'MarkerSize', 8);
xlabel('Mesh Size (N = M)');
ylabel('Number of Iterations');
title('Solver Convergence Iterations vs Mesh Refinement (Case 1)');
legend('Case 1');
grid on;
set(gca, 'XScale', 'log');

% Solver Residual vs Mesh Size
figure('Name', 'Solver Residual vs Mesh Refinement (Case 1)', 'NumberTitle', 'off', 'Position', figPos);
loglog(gridSizes, convergenceMetrics.solver_residual(1, :), ...
       'd-', 'LineWidth', 2.5, 'MarkerSize', 8);
xlabel('Mesh Size (N = M)');
ylabel('Final Residual');
title('Solver Residual vs Mesh Refinement (Case 1)');
legend('Case 1');
grid on;
set(gca, 'XScale', 'log');

%% Summary Table
fprintf('\n=== CONVERGENCE SUMMARY TABLE ===\n');
fprintf('Grid  | Spacing  | L2-Error   | L∞-Error   | Iterations | Residual\n');
fprintf('------|----------|------------|------------|------------|----------\n');
for gridIdx = 1:numGrids
    fprintf('%3d   | %.2e | %.2e | %.2e |    %d      | %.2e\n', ...
        gridSizes(gridIdx), convergenceMetrics.spacing_dx(gridIdx), ...
        convergenceMetrics.error_L2(1, gridIdx), ...
        convergenceMetrics.error_Linf(1, gridIdx), ...
        convergenceMetrics.solver_iterations(1, gridIdx), ...
        convergenceMetrics.solver_residual(1, gridIdx));
end

%% Save Results
save('convergence_analysis_results.mat', 'convergenceMetrics');
fprintf('\n\nAnalysis results saved to: convergence_analysis_results.mat\n');
fprintf('Convergence validation study completed successfully!\n');