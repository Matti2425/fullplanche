function [C_D, C_L, Gamma, D, L, p_cyl] = final_calc(x_P, y_P, p_P, v_xP, v_yP, x_cyl, y_cyl, theta, Geometry, w, v_in)
% FINAL_CALC Calculates aerodynamic properties of a cylinder in flow
%
% This function calculates the aerodynamic properties of a rotating cylinder immersed in 
% a potential flow. It determines the lift and drag forces, their respective coefficients, 
% and the circulation around the cylinder.
%
% INPUTS:
%   x_P, y_P  - 1D arrays of grid point coordinates
%   p_P       - Pressure field at each computational node (M+2, N+2)
%   v_xP, v_yP - Velocity components at each grid point (M+2, N+2)
%   x_cyl, y_cyl - Cylinder surface coordinates
%   theta     - Angular positions along the cylinder surface
%   Geometry  - Structure containing ObjectCenter, D (diameter), rho_ref
%   w         - Angular velocity of the cylinder (rad/s)
%   v_in      - Free-stream velocity
%
% OUTPUTS:
%   C_D   - Drag coefficient
%   C_L   - Lift coefficient
%   Gamma - Circulation around the cylinder
%   D     - Drag force per unit depth
%   L     - Lift force per unit depth
%   p_cyl - Interpolated pressure on cylinder surface
%--------------------------------------------------------------------------

    % Extract geometry parameters
    R_cyl = Geometry.D / 2;  % Radius
    rho_ref = Geometry.rho_ref;
    
    % Interpolate pressure on the cylinder surface
    p_cyl = interp2(x_P, y_P, p_P, x_cyl, y_cyl, 'spline');
    
    % Interpolate velocity components at the cylinder surface
    v_x_cyl = interp2(x_P, y_P, v_xP, x_cyl, y_cyl, 'spline');
    v_y_cyl = interp2(x_P, y_P, v_yP, x_cyl, y_cyl, 'spline');
    
    % Compute tangential velocity (vectorized)
    v_theta = v_x_cyl .* cos(theta) + v_y_cyl .* sin(theta);
    
    % Compute circulation
    Gamma = 2 * pi * R_cyl^2 * w;
    
    % Compute lift and drag forces (vectorized)
    L = -6.5 * trapz(theta, p_cyl .* sin(theta)) * R_cyl;
    D = -trapz(theta, p_cyl .* cos(theta)) * R_cyl;
    
    % Compute lift and drag coefficients
    C_D = D / (0.5 * rho_ref * v_in^2 * (2 * R_cyl));
    C_L = L / (0.5 * rho_ref * v_in^2 * (2 * R_cyl));
    
    % Plot pressure distribution
    figure;
    plot(theta, p_cyl);
    xlabel('\theta'); 
    ylabel('Interpolated Pressure');
    title('Pressure over cylinder surface');
    
    % Display results
    fprintf('Lift Coefficient (C_L): %f\n', C_L);
    fprintf('Drag Coefficient (C_D): %f\n', C_D);
    fprintf('Circulation (Gamma): %f\n', Gamma);
    
end