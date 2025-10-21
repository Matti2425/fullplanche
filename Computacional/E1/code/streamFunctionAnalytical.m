function psiAnalytical = streamFunctionAnalytical(gridX, gridY, domain, velocity, simCase, rotSpeed, obstacleRadius)
% STREAMFUNCTIONANALYTICAL Computes analytical stream function solutions
%
% Calculates the analytical stream function for three different flow configurations:
% - Case 1: Uniform channel flow
% - Case 2: Static obstacle with potential flow
% - Case 3: Rotating obstacle with circulation
%
% INPUTS:
%   gridX, gridY       - 1D mesh coordinate arrays
%   domain             - Geometry structure with center coordinates and dimensions
%   velocity           - Inlet flow velocity
%   simCase            - Case selector (1, 2, or 3)
%   rotSpeed           - Angular velocity for rotating case
%   obstacleRadius     - Radius of cylindrical obstacle
%
% OUTPUTS:
%   psiAnalytical - Analytical stream function field (ny+2, nx+2)
%--------------------------------------------------------------------------

    ny = length(gridY) - 2;
    nx = length(gridX) - 2;
    psiAnalytical = zeros(ny+2, nx+2);

    centerX = domain.center_x;
    centerY = domain.center_y;
    domainHeight = domain.height;
    domainLength = domain.length;

    % Create mesh grid for computation
    [meshX, meshY] = meshgrid(gridX, gridY);

    switch simCase
        case 1
            % CASE 1: Uniform channel flow (simple linear profile)
            % psi(x, y) = velocity * y
            psiAnalytical = velocity * meshY;

        case 2
            % CASE 2: Potential flow around static obstacle (uniform + doublet)
            % psi = psi_uniform + psi_doublet
            
            % Distance from obstacle center
            radiusSq = (meshX - centerX).^2 + (meshY - centerY).^2;
            radius = sqrt(radiusSq);

            % Uniform flow component
            psiUniform = velocity * meshY;

            % Doublet strength representing the obstacle effect
            doubletStrength = 2 * pi * obstacleRadius^2 * velocity;

            % Doublet component (only outside obstacle)
            psiDoublet = zeros(size(meshX));
            outerRegion = radius > obstacleRadius;
            psiDoublet(outerRegion) = (doubletStrength / (2*pi)) * ...
                                     (meshY(outerRegion) - centerY) ./ radiusSq(outerRegion);

            % Total solution: uniform + doublet
            psiAnalytical = psiUniform + psiDoublet;

            % Inside obstacle: constant stream function
            psiAnalytical(radius <= obstacleRadius) = velocity * centerY;

        case 3
            % CASE 3: Rotating obstacle with circulation (static + vortex)
            % psi = psi_static + psi_circulation

            radius = sqrt((meshX - centerX).^2 + (meshY - centerY).^2);
            radiusSq = radius.^2;

            % Static component (uniform + doublet)
            doubletStrength = 2 * pi * obstacleRadius^2 * velocity;
            psiUniform = velocity * meshY;

            psiStatic = psiUniform;
            outerRegion = radius > obstacleRadius;
            psiStatic(outerRegion) = psiStatic(outerRegion) + ...
                                    (doubletStrength / (2*pi)) * ...
                                    (meshY(outerRegion) - centerY) ./ radiusSq(outerRegion);
            psiStatic(radius <= obstacleRadius) = velocity * centerY;

            % Circulation component from rotation
            circulationStrength = 2 * pi * obstacleRadius^2 * rotSpeed;
            psiCirculation = zeros(size(meshX));
            psiCirculation(outerRegion) = -(circulationStrength / (2*pi)) * ...
                                         log(radius(outerRegion) / obstacleRadius);

            % Total solution: static + circulation effect
            psiAnalytical = psiStatic + psiCirculation;
    end

end