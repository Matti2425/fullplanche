%% CFD Simulation - Potential Flow Around a Cylinder
% E1 Implementation
clearvars; close all; clc;

%% Read Input Data
[Geometry, v_in, p_in, T_in, psi_body, gas, Numerical, Init, w, CASE] = readInputData();

%% Generate Mesh
[x_P, y_P] = generateMesh(Geometry, Numerical);

%% Generate Cylinder Geometry and Solid Matrix
[x_cyl, y_cyl, S, theta] = solidMaterial(Geometry, x_P, y_P, Numerical);
if CASE == 1
    S = S*0;
end
%% Initialize Fields
[psi, rho] = initParams(Init, Numerical, S);

%% Evaluate Coefficients
[a_E, a_W, a_N, a_S, b_P, a_P] = evalCoefficients(psi, psi_body, Numerical, Geometry, v_in, x_P, y_P, S, Init);

%% Solve System (Gauss-Seidel Iteration)
fprintf('Starting Gauss-Seidel solver...\n');
[psi, i, error] = solveSystem(psi, psi_body, S, a_E, a_W, a_N, a_S, a_P, b_P, v_in, y_P, Numerical, Geometry);
fprintf('Converged after %d iterations with error = %.2e\n', i, error);

%% Compute Velocities
[v_x, v_y, v_mag] = compute_vel(psi, x_P, y_P, Numerical, v_in, Geometry.rho_ref, rho);

%% Compute Thermodynamic Properties
[T, P, rho, cp_mean] = compute_thermo(v_mag, Numerical, Init, gas);

%% Compute Aerodynamic Coefficients
if CASE ~= 1
    [C_D, C_L, Gamma, D, L, p_cyl] = final_calc(x_P, y_P, P, v_x, v_y, x_cyl, y_cyl, theta, Geometry, w, v_in);
    
    % Display results
    fprintf('\n=== Aerodynamic Results ===\n');
    fprintf('Drag Coefficient (C_D): %.6f\n', C_D);
    fprintf('Lift Coefficient (C_L): %.6f\n', C_L);
    fprintf('Circulation (Gamma): %.6f\n', Gamma);
    fprintf('Drag Force (D): %.6f N\n', D);
    fprintf('Lift Force (L): %.6f N\n', L);
end

%% Plotting
% Calculate grid spacing
dx = Geometry.L / Numerical.N;
dy = Geometry.H / Numerical.M;

figure('Name', 'Stream Function Contours', 'NumberTitle', 'off');
plotStreamLines(psi, x_P, y_P);
if CASE ~= 1
    hold on;
    plot(x_cyl, y_cyl, 'r-', 'LineWidth', 2);
    hold off;
end
xlabel('x [m]'); ylabel('y [m]');
title('Stream Function Contours');
axis equal tight;
name = "figures/" + "Case_"+ CASE + "/" + CASE + "_Streamfunc.png";
exportgraphics(gca, name, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'none');


figure('Name', 'Pressure Field (Pa)', 'NumberTitle', 'off');
plotCFD('Pressure (Pa)', P, Numerical.N, Numerical.M, dx, dy);
if CASE ~= 1
    hold on;
    plot(x_cyl, y_cyl, 'k-', 'LineWidth', 2);
    hold off;
end
name = "figures/" + "Case_"+ CASE + "/" + CASE + "_Pressure.png";
exportgraphics(gca, name, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'none');

figure('Name', 'Velocity Magnitude (m/s)', 'NumberTitle', 'off');
plotCFD('Velocity (m/s)', v_mag, Numerical.N, Numerical.M, dx, dy);
if CASE ~= 1
    hold on;
    plot(x_cyl, y_cyl, 'k-', 'LineWidth', 2);
    hold off;
end
name = "figures/" + "Case_"+ CASE + "/" + CASE + "_Velocity.png";
exportgraphics(gca, name, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'none');

if CASE ~= 1
    figure('Name', 'Pressure Distribution on Cylinder', 'NumberTitle', 'off');
    plot(theta * 180/pi, p_cyl, 'b-', 'LineWidth', 2);
    xlabel('Angle (degrees)');
    ylabel('Pressure (Pa)');
    title('Pressure Distribution on Cylinder Surface');
    grid on;
end

fprintf('\nSimulation completed successfully!\n');