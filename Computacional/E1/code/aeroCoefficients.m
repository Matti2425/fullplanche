function [coeffDrag, coeffLift, circulation, forceDrag, forceLift, pressureOnBody] = ...
         aeroCoefficients(gridX, gridY, pressField, velX, velY, cylX, cylY, angularPos, domain, velocity)
% AEROCOEFFICIENTS Computes aerodynamic performance metrics
%
% Calculates drag and lift coefficients, forces, and circulation
% around a cylindrical obstacle using pressure and velocity fields.
% Uses Kutta-Joukowski theorem: L = rho * V_inf * Gamma
%
% INPUTS:
%   gridX, gridY     - Coordinate vectors
%   pressField       - Pressure field at grid points
%   velX, velY       - Velocity components at grid points
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

    % Compute tangential velocity at surface (counter-clockwise positive)
    % Normal vector points outward from center: (cos(theta), sin(theta))
    % Tangent vector (counter-clockwise): (-sin(theta), cos(theta))
    velTangential = -velX_surface .* sin(angularPos) + velY_surface .* cos(angularPos);

    % Compute circulation from velocity line integral
    % Gamma = ∮ V·dl around cylinder
    circulation = trapz(angularPos, velTangential) * obstacleRadius;

    % Compute aerodynamic forces using pressure integration
    % Normal vector: n = (cos(theta), sin(theta))
    % Pressure force: F = -∫ p * n dS
    % Components: Fx = -∫ p*cos(theta) dS, Fy = -∫ p*sin(theta) dS
    
    forceDrag = -trapz(angularPos, pressureOnBody .* cos(angularPos)) * obstacleRadius;
    forceLift = -trapz(angularPos, pressureOnBody .* sin(angularPos)) * obstacleRadius;

    % Reference dynamic pressure for non-dimensionalization
    qRef = 0.5 * rho_reference * velocity^2;
    referenceLength = domain.diameter;  % Diameter is characteristic length

    % Compute dimensionless coefficients
    coeffDrag = forceDrag / (qRef * referenceLength);
    coeffLift = forceLift / (qRef * referenceLength);

end