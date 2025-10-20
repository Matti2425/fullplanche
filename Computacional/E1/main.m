%% Potential Flow Simulation Around Cylindrical Obstacle
% Gerard Martínez Vilella
%
clearvars; close all; clc;

%% Select case
simCase = 2; % 1- Uniform flow 2- Static Cylinder 3- Rotating Cylinder

%% Configuration and Parameter Setup
[domain, inlet, boundary, material, solver, state, rotSpeed] = setupSimulation(simCase);

%% Mesh Generation
[gridX, gridY] = createComputationalGrid(domain, solver);

%% Obstacle Geometry Definition
[cylX, cylY, fluidMask, angularPos] = defineObstacles(domain, gridX, gridY, solver, simCase);

%% Field Initialization
[streamFunc, density] = initializeFields(state, solver, fluidMask);

%% Construct Discretization Matrix
[coeff_east, coeff_west, coeff_north, coeff_south, coeff_center, rhs] = ...
    constructMatrix(streamFunc, boundary.psi_surface, solver, domain, inlet.velocity, ...
                    gridX, gridY, fluidMask, state);

%% Iterative Solution (Gauss-Seidel Method)
fprintf('Initiating iterative solver...\n');
[streamFunc, numIter, residual] = iterativeSolver(streamFunc, boundary.psi_surface, fluidMask, ...
                                                   coeff_east, coeff_west, coeff_north, coeff_south, ...
                                                   coeff_center, rhs, inlet.velocity, gridY, solver, domain);
fprintf('Convergence achieved in %d iterations with residual = %.2e\n', numIter, residual);

%% Velocity Field Reconstruction
[velX, velY, velMag] = velocityField(streamFunc, gridX, gridY, solver, inlet.velocity, ...
                                     domain.ref_density, density);

%% Thermodynamic Property Calculation
[tempField, pressField, rhoField, cpAvg] = thermoProperties(velMag, solver, state, material);

%% Aerodynamic Coefficient Evaluation
if simCase ~= 1
    [coeffDrag, coeffLift, circulation, forceDrag, forceLift, pressureOnBody] = ...
        aeroCoefficients(gridX, gridY, pressField, velX, velY, cylX, cylY, ...
                        angularPos, domain, rotSpeed, inlet.velocity);
    
    % Results Display
    fprintf('\n=== Aerodynamic Performance Metrics ===\n');
    fprintf('Drag Coefficient (C_D):  %.6f\n', coeffDrag);
    fprintf('Lift Coefficient (C_L):  %.6f\n', coeffLift);
    fprintf('Circulation (Gamma):     %.6f\n', circulation);
    fprintf('Drag Force (D):          %.6f N\n', forceDrag);
    fprintf('Lift Force (L):          %.6f N\n', forceLift);
end

%% Visualization
% Compute mesh spacing for proper scaling
dx = domain.length / solver.nx;
dy = domain.height / solver.ny;

figure('Name', 'Flow Pattern Analysis', 'NumberTitle', 'off');
displayFlowPattern(streamFunc, gridX, gridY);
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

fprintf('\nSimulation execution completed successfully!\n');