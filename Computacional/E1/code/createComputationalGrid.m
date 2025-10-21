function [gridX, gridY] = createComputationalGrid(domain, solver)
% CREATECOMPUTATIONALGRID Generates a structured Cartesian computational mesh
%
% Constructs a uniform grid for a 2D domain with specified dimensions.
% Grid points are positioned at cell centers and boundaries.
%
% Inputs:
%   domain - Structure containing length and height of the domain
%   solver - Structure containing nx and ny (number of grid points)
%
% Outputs:
%   gridX - X-coordinate vector (includes boundaries and interior points)
%   gridY - Y-coordinate vector (includes boundaries and interior points)
%--------------------------------------------------------------------------

    nx = solver.nx;
    ny = solver.ny;
    domainLength = domain.length;
    domainHeight = domain.height;

    % Generate X-coordinate array
    gridX = zeros(1, nx+2);
    gridX(2:end-1) = (0.5:1:nx-0.5) * (domainLength/nx);  % Interior cell centers
    gridX(1) = 0;                                           % Left boundary
    gridX(end) = domainLength;                              % Right boundary

    % Generate Y-coordinate array
    gridY = zeros(1, ny+2);
    gridY(2:end-1) = (0.5:1:ny-0.5) * (domainHeight/ny);  % Interior cell centers
    gridY(1) = 0;                                           % Bottom boundary
    gridY(end) = domainHeight;                              % Top boundary

end