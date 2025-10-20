function [coeffDrag, coeffLift, circulation, forceDrag, forceLift, pressureOnBody] = ...
         aeroCoefficients(gridX, gridY, pressField, velX, velY, cylX, cylY, angularPos, domain, rotSpeed, velocity)
% AEROCOEFFICIENTS Computes aerodynamic performance metrics
%
% Calculates drag and lift coefficients, forces, and circulation
% around a cylindrical obstacle using pressure and velocity fields.
%
% INPUTS:
%   gridX, gridY     - Coordinate vectors
%   pressField       - Pressure field at grid points (ny+2, nx+2)
%   velX, velY       - Velocity components at grid points (ny+2, nx+2)
%   cylX, cylY       - Cylinder surface coordinate arrays
%   angularPos       - Angular positions around cylinder
%   domain           - Geometry structure (center, diameter, rho_ref)
%   rotSpeed         - Angular velocity (rad/s)
%   velocity         - Free-stream velocity
%
% OUTPUTS:
%   coeffDrag   - Drag coefficient
%   coeffLift   - Lift coefficient
%   circulation - Circulation around cylinder (Gamma)
%   forceDrag   - Drag force per unit depth
%   forceLift   - Lift force per unit depth
%   pressureOnBody - Pressure distribution on surface
%--------------------------------------------------------------------------

    % Extract obstacle parameters
    obstacleRadius = domain.diameter / 2;
    rho_reference = domain.ref_density;

    % Interpolate pressure on cylinder surface
    pressureOnBody = interp2(gridX, gridY, pressField, cylX, cylY, 'spline');

    % Interpolate velocity components on cylinder surface
    velX_surface = interp2(gridX, gridY, velX, cylX, cylY, 'spline');
    velY_surface = interp2(gridX, gridY, velY, cylX, cylY, 'spline');

    % Compute tangential velocity at surface
    velTangential = velX_surface .* cos(angularPos) + velY_surface .* sin(angularPos);

    % Compute circulation using Stokes theorem
    circulation = 2 * pi * obstacleRadius^2 * rotSpeed;

    % Compute aerodynamic forces using pressure integration
    % Forces per unit depth (span-wise dimension = 1)
    forceLift = -trapz(angularPos, pressureOnBody .* sin(angularPos)) * obstacleRadius;
    forceDrag = -trapz(angularPos, pressureOnBody .* cos(angularPos)) * obstacleRadius;

    % Reference dynamic pressure for non-dimensionalization
    qRef = 0.5 * rho_reference * velocity^2;
    referenceLength = 2 * obstacleRadius;  % Diameter is characteristic length

    % Compute dimensionless coefficients
    coeffDrag = forceDrag / (qRef * referenceLength);
    coeffLift = forceLift / (qRef * referenceLength);

end