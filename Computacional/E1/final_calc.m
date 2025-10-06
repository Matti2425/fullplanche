function [C_D, C_L, Gamma, D, L, p_cyl, theta] = final_calc(X, Y, P, v_x, v_y, Geometry, Vin, gas)
% FINAL_CALC Calculates aerodynamic properties of a cylinder in flow
%
% This function calculates the aerodynamic properties of a cylinder immersed in 
% a flow field. It determines the lift and drag forces, their respective coefficients, 
% and the circulation around the cylinder.
%
% INPUTS:
%   X, Y     - Grid point coordinates in the computational domain (N+2, M+2)
%   P        - Pressure field at each computational node (N+2, M+2)
%   v_x, v_y - Velocity components at each grid point (N+2, M+2)
%   Geometry - Structure containing:
%              * ObjectCenter - [x_center, y_center] of cylinder
%              * D - Diameter of cylinder
%              * rho_ref - Reference density
%   Vin      - Free-stream velocity
%   gas      - Structure containing gas properties (optional, for future use)
%
% OUTPUTS:
%   C_D    - Drag coefficient
%   C_L    - Lift coefficient
%   Gamma  - Circulation around the cylinder
%   D      - Drag force per unit depth
%   L      - Lift force per unit depth
%   p_cyl  - Interpolated pressure on cylinder surface
%   theta  - Angular positions along the cylinder surface
%--------------------------------------------------------------------------

    % Extract geometry parameters
    x_center = Geometry.ObjectCenter(1);
    y_center = Geometry.ObjectCenter(2);
    R_cyl = Geometry.D / 2;  % Radius
    rho_ref = Geometry.rho_ref;
    
    % Generate cylinder surface coordinates
    n_points = 200;  % Number of points around cylinder
    theta = linspace(0, 2*pi, n_points);
    
    % Cylinder surface coordinates (vectorized)
    % Sample slightly outside the cylinder to avoid solid region
    R_sample = R_cyl * 1.05;  % 5% larger radius to sample in fluid region
    x_cyl = x_center + R_sample .* cos(theta);
    y_cyl = y_center + R_sample .* sin(theta);
    
    % Remove NaN values from the data before interpolation
    % (NaN can occur in solid regions)
    valid_idx = ~isnan(P(:)) & ~isnan(v_x(:)) & ~isnan(v_y(:));
    
    % Create scattered interpolants for better performance and robustness
    F_p = scatteredInterpolant(X(valid_idx), Y(valid_idx), P(valid_idx), 'natural', 'nearest');
    F_vx = scatteredInterpolant(X(valid_idx), Y(valid_idx), v_x(valid_idx), 'natural', 'nearest');
    F_vy = scatteredInterpolant(X(valid_idx), Y(valid_idx), v_y(valid_idx), 'natural', 'nearest');
    
    % Interpolate pressure on the cylinder surface
    p_cyl = F_p(x_cyl, y_cyl);
    
    % Interpolate velocity components at the cylinder surface
    v_x_cyl = F_vx(x_cyl, y_cyl);
    v_y_cyl = F_vy(x_cyl, y_cyl);
    
    % Compute tangential velocity (vectorized)
    % v_theta = v_x * cos(theta) + v_y * sin(theta)
    v_theta = v_x_cyl .* cos(theta) + v_y_cyl .* sin(theta);
    
    % Compute circulation using line integral of velocity
    % Gamma = integral of v_theta * R * dtheta around the cylinder
    Gamma = trapz(theta, v_theta) * R_cyl;
    
    % Alternative: For rotating cylinder, Gamma = 2*pi*R^2*omega
    % If you have angular velocity (omega), uncomment:
    % if isfield(Geometry, 'omega')
    %     Gamma = 2 * pi * R_cyl^2 * Geometry.omega;
    % end
    
    % Compute lift and drag forces per unit depth (vectorized)
    % Force components from pressure integration:
    % F_x = -integral(p * cos(theta) * R * dtheta)  [Drag]
    % F_y = -integral(p * sin(theta) * R * dtheta)  [Lift]
    
    D = -trapz(theta, p_cyl .* cos(theta)) * R_cyl;
    L = -trapz(theta, p_cyl .* sin(theta)) * R_cyl;
    
    % Compute dynamic pressure for normalization
    q_inf = 0.5 * rho_ref * Vin^2;
    
    % Compute lift and drag coefficients
    % Reference area for 2D cylinder = diameter * unit depth = D * 1
    A_ref = 2 * R_cyl;  % Diameter (per unit depth)
    
    C_D = D / (q_inf * A_ref);
    C_L = L / (q_inf * A_ref);
    
    % Display results
    fprintf('\n========== AERODYNAMIC RESULTS ==========\n');
    fprintf('Lift Force (L):          %12.6f N/m\n', L);
    fprintf('Drag Force (D):          %12.6f N/m\n', D);
    fprintf('Lift Coefficient (C_L):  %12.6f\n', C_L);
    fprintf('Drag Coefficient (C_D):  %12.6f\n', C_D);
    fprintf('Circulation (Gamma):     %12.6f m²/s\n', Gamma);
    fprintf('=========================================\n\n');
    
    % Plot pressure distribution around cylinder
    figure('Name', 'Pressure Distribution on Cylinder Surface');
    plot(theta * 180/pi, p_cyl, 'b-', 'LineWidth', 2);
    xlabel('Angular Position \theta (degrees)', 'FontSize', 12);
    ylabel('Pressure (Pa)', 'FontSize', 12);
    title('Pressure Distribution on Cylinder Surface', 'FontSize', 14);
    grid on;
    xlim([0 360]);
    
    % Add reference lines
    hold on;
    yline(mean(p_cyl), 'r--', 'Mean Pressure', 'LineWidth', 1.5);
    hold off;
    
    % Optional: Plot pressure coefficient
    figure('Name', 'Pressure Coefficient Distribution');
    C_p = (p_cyl - mean(p_cyl(:))) / q_inf;
    plot(theta * 180/pi, C_p, 'r-', 'LineWidth', 2);
    xlabel('Angular Position \theta (degrees)', 'FontSize', 12);
    ylabel('Pressure Coefficient C_p', 'FontSize', 12);
    title('Pressure Coefficient Distribution on Cylinder', 'FontSize', 14);
    grid on;
    xlim([0 360]);
    yline(0, 'k--', 'LineWidth', 1);
    
end