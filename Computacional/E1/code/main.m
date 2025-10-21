%% Potential Flow Simulation Around Cylindrical Obstacle
% Gerard Martínez Vilella
%
clearvars; close all;

%% Select case
simCase = 3; % 1- Uniform flow 2- Static Cylinder 3- Rotating Cylinder
PLOTS = true;

%% Configuration and Parameter Setup
[domain, inlet, boundary, material, solver, state, rotSpeed] = setupSimulation(simCase);

%% Mesh Generation
[gridX, gridY] = createComputationalGrid(domain, solver);

%% Obstacle Geometry Definition
[cylX, cylY, fluidMask, angularPos] = defineObstacles(domain, gridX, gridY, solver, simCase);

%% Field Initialization
[psi, density] = initializeFields(state, solver, fluidMask);

%% Construct Discretization Matrix
[a_E, a_W, a_N, a_S, a_P, rhs] = ...
    constructMatrix(psi, boundary.psi_surface, solver, domain, inlet.velocity, ...
                    gridX, gridY, fluidMask, state);

%% Iterative Solution (Gauss-Seidel Method)
fprintf('Initiating iterative solver...\n');
[psi, numIter, residual] = iterativeSolver(psi, boundary.psi_surface, fluidMask, ...
                                                   a_E, a_W, a_N, a_S, ...
                                                   a_P, rhs, inlet.velocity, gridY, solver, domain);
fprintf('Convergence achieved in %d iterations with residual = %.2e\n', numIter, residual);

%% Velocity Field Reconstruction
[velX, velY, velMag] = velocityField(psi, gridX, gridY, solver, inlet.velocity, ...
                                     domain.ref_density, density);
%% Verification for Case 1
if simCase == 1
    % Check if psi is uniform in x-direction and linear in y
    fprintf('\n=== Case 1 Verification: Stream Function ===\n');

    % Sample psi at middle y-location
    y_mid = ceil((solver.ny+2)/2);
    psi_x_profile = psi(y_mid, :);
    psi_variation_x = max(psi_x_profile) - min(psi_x_profile);
    fprintf('psi variation in x-direction: %.2e \n', psi_variation_x);

    % Check linearity with y at inlet
    psi_at_inlet = psi(:, 50);
    y_coords = linspace(0, domain.height, solver.ny+2);
    psi_analytical = inlet.velocity * y_coords;
    error_psi = max(abs(psi_at_inlet - psi_analytical'));
    fprintf('Max error from analytical psi(y): %.2e \n', error_psi);

    % Display some sample values
    fprintf('\nSample psi values at inlet (i=1):\n');
    for j = 1:30:solver.ny+2
        fprintf('  y=%.3f: psi_numerical=%.6f, psi_analytical=%.6f\n', ...
            y_coords(j), psi_at_inlet(j), psi_analytical(j));
    end
end

%% Thermodynamic Property Calculation
[tempField, pressField, rhoField, cpAvg] = thermoProperties(velMag, solver, state, material);

%% Aerodynamic Coefficient Evaluation
if simCase ~= 1
    [coeffDrag, coeffLift, circulation, forceDrag, forceLift, pressureOnBody] = aeroCoefficients(gridX, gridY, pressField, velX, velY, cylX, cylY, ...
                                                                                                    angularPos, domain, inlet.velocity);
    
    % Results Display
    fprintf('\n=== Aerodynamic Performance Metrics ===\n');
    fprintf('Drag Coefficient (C_D):  %.6f\n', coeffDrag);
    fprintf('Lift Coefficient (C_L):  %.6f\n', coeffLift);
    fprintf('Circulation (Gamma):     %.6f\n', circulation);
    fprintf('Drag Force (D):          %.6f N\n', forceDrag);
    fprintf('Lift Force (L):          %.6f N\n', forceLift);
end

if PLOTS
    %% Visualization
    % Compute mesh spacing for proper scaling
    dx = domain.length / solver.nx;
    dy = domain.height / solver.ny;
    
    figure('Name', 'Flow Pattern Analysis', 'NumberTitle', 'off');
    displayFlowPattern(psi, gridX, gridY);
    if simCase ~= 1
        hold on;
        plot(cylX, cylY, 'r-', 'LineWidth', 2);
        hold off;
    end
    xlabel('x [m]'); ylabel('y [m]');
    title('Stream Function Streamlines');
    axis equal tight;
    name = "figures/" + "Case_"+ simCase + "/" + simCase + "_Streamfunc.png";
    exportgraphics(gca, name, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'none');
    
    
    figure('Name', 'Pressure Distribution (Pa)', 'NumberTitle', 'off');
    visualizeScalarField('Pressure (Pa)', pressField, solver.nx, solver.ny, dx, dy);
    if simCase ~= 1
        hold on;
        plot(cylX, cylY, 'k-', 'LineWidth', 2);
        hold off;
    end
    name = "figures/" + "Case_"+ simCase + "/" + simCase + "_Pressure.png";
    exportgraphics(gca, name, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'none');
    
    figure('Name', 'Temperature Profile (K)', 'NumberTitle', 'off');
    visualizeScalarField('Temperature (K)', tempField, solver.nx, solver.ny, dx, dy);
    if simCase ~= 1
        hold on;
        plot(cylX, cylY, 'k-', 'LineWidth', 2);
        hold off;
    end
    name = "figures/" + "Case_"+ simCase + "/" + simCase + "_Temperature.png";
    exportgraphics(gca, name, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'none');
    
    figure('Name', 'Velocity Field (m/s)', 'NumberTitle', 'off');
    visualizeScalarField('Velocity (m/s)', velMag, solver.nx, solver.ny, dx, dy);
    if simCase ~= 1
        hold on;
        plot(cylX, cylY, 'k-', 'LineWidth', 2);
        hold off;
    end
    name = "figures/" + "Case_"+ simCase + "/" + simCase + "_Velocity.png";
    exportgraphics(gca, name, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'none');
    
    if simCase ~= 1
        figure('Name', 'Surface Pressure Distribution', 'NumberTitle', 'off');
        plot(angularPos * 180/pi, pressureOnBody, 'b-', 'LineWidth', 2);
        xlabel('Angular Position (degrees)');
        ylabel('Pressure (Pa)');
        title('Pressure Variation Around Obstacle Surface');
        grid on;
    end
end

fprintf('\nSimulation execution completed successfully!\n');