function [cylX, cylY, fluidMask, angularPos] = defineObstacles(domain, gridX, gridY, solver, simCase)
% DEFINEOBSTACLES Defines cylindrical obstacle geometry in the domain
%
% Creates a circular obstacle at domain center and marks fluid/solid regions.
%
% Inputs:
%   domain    - Structure with center coordinates and diameter
%   gridX, gridY - Grid coordinate vectors
%   solver    - Structure containing nx and ny
%
% Outputs:
%   cylX, cylY   - Coordinates defining obstacle perimeter
%   fluidMask    - Binary matrix (1=fluid, 0=solid)
%   angularPos   - Angular positions around obstacle surface
%--------------------------------------------------------------------------

    nx = solver.nx;
    ny = solver.ny;
    obstacleX = domain.center_x;
    obstacleY = domain.center_y;
    obstacleRadius = domain.diameter / 2;

    % Initialize fluid region marker (1 = fluid everywhere initially)
    fluidMask = ones(ny+2, nx+2);

    % Generate obstacle perimeter points
    angularPos = linspace(0, 2*pi, 10000);
    cylX = obstacleX + obstacleRadius * cos(angularPos);
    cylY = obstacleY + obstacleRadius * sin(angularPos);

    % Mark solid region (points inside obstacle)
    if simCase ~= 1
        for j = 1:ny+2
            for i = 1:nx+2
                distanceToCenter = (gridX(i) - obstacleX)^2 + (gridY(j) - obstacleY)^2;
                if distanceToCenter <= obstacleRadius^2
                    fluidMask(j, i) = 0;  % Mark as solid
                end
            end
        end
    end

end