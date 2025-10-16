function [nod_x, nod_y, center_x, center_y, sx_w, sx_e, sx_p, sy_s, sy_n, sy_p] = generateMesh(Geometry, Numerical)
% GENERATEMESH Creates a structured mesh for convection-diffusion problem
%
% This function generates a uniform computational grid with cell centers,
% face positions, and geometric parameters needed for finite volume discretization.
%
% INPUTS:
%   Geometry  - Structure containing L, H (domain dimensions)
%   Numerical - Structure containing N, M (number of cells)
%
% OUTPUTS:
%   nod_x, nod_y       - Node coordinates (cell centers) including boundaries
%   center_x, center_y - Cell face centers
%   sx_w, sx_e, sx_p   - Geometric parameters in x-direction (west, east, P)
%   sy_s, sy_n, sy_p   - Geometric parameters in y-direction (south, north, P)
%--------------------------------------------------------------------------

    L = Geometry.L;
    H = Geometry.H;
    N = Numerical.N;
    M = Numerical.M;
    
    % Cell dimensions
    A_x = L / N;
    A_y = H / M;
    
    % Generate cell face centers (0:N and 0:M)
    center_x = (0:N) * A_x;
    center_y = (0:M) * A_y;
    
    % Generate cell centers (interior nodes)
    nodvec_x = ((0:N-1) * A_x) + (A_x / 2);
    nodvec_y = ((0:M-1) * A_y) + (A_y / 2);
    
    % Add boundary nodes (ghost cells)
    nod_x = [0, nodvec_x, L];
    nod_y = [0, nodvec_y, H];
    
    % Preallocate geometric parameter arrays
    sx_w = zeros(size(nod_x));
    sx_e = zeros(size(nod_x));
    sx_p = zeros(size(nod_x));
    sy_n = zeros(size(nod_y));
    sy_s = zeros(size(nod_y));
    sy_p = zeros(size(nod_y));
    
    % X-direction geometric parameters (vectorized)
    % Interior cells
    sx_w(2:end-1) = nod_x(2:end-1) - nod_x(1:end-2);
    sx_e(2:end-1) = nod_x(3:end) - nod_x(2:end-1);
    
    % Ensure center_x matches nod_x size
    if numel(center_x) ~= numel(nod_x)
        center_x = linspace(min(nod_x), max(nod_x), numel(nod_x));
    end
    sx_p(2:end-1) = center_x(2:end-1) - center_x(1:end-2);
    
    % X-direction boundary conditions
    sx_w(1) = 0;
    sx_w(end) = nod_x(end) - nod_x(end-1);
    sx_e(1) = nod_x(2) - nod_x(1);
    sx_e(end) = 0;
    sx_p(1) = center_x(2) / 2;
    sx_p(end) = (center_x(end) - center_x(end-1)) / 2;
    
    % Y-direction geometric parameters (vectorized)
    % Interior cells
    sy_s(2:end-1) = nod_y(2:end-1) - nod_y(1:end-2);
    sy_n(2:end-1) = nod_y(3:end) - nod_y(2:end-1);
    
    % Ensure center_y matches nod_y size
    if numel(center_y) ~= numel(nod_y)
        center_y = linspace(min(nod_y), max(nod_y), numel(nod_y));
    end
    sy_p(2:end-1) = center_y(2:end-1) - center_y(1:end-2);
    
    % Y-direction boundary conditions
    sy_s(1) = 0;
    sy_s(end) = nod_y(end) - nod_y(end-1);
    sy_n(1) = nod_y(2) - nod_y(1);
    sy_n(end) = 0;
    sy_p(1) = center_y(2) / 2;
    sy_p(end) = (center_y(end) - center_y(end-1)) / 2;

end