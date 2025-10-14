function [x_P, y_P] = generateMesh(Geometry, NumericalParameters)
% MESH_GENERATOR Creates a structured Cartesian mesh for a 2D domain.
%
% This function generates a uniform computational grid for a domain 
% of dimensions L × H, discretized with N points in the x-direction 
% and M points in the y-direction.
%
% Inputs:
%   Geometry            - Structure containing L (length) and H (height)
%   NumericalParameters - Structure containing N and M (grid points)
%
% Outputs:
%   x_P - 1D array containing the x-coordinates of the cell centers and boundaries
%   y_P - 1D array containing the y-coordinates of the cell centers and boundaries
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    L = Geometry.L;
    H = Geometry.H;

    % Initialize 1D coordinate arrays
    x_P = zeros(1, N+2);
    x_P(2:end-1) = (0.5:1:N-0.5) * (L/N);
    x_P(1) = 0;
    x_P(end) = L;

    y_P = zeros(1, M+2);
    y_P(2:end-1) = (0.5:1:M-0.5) * (H/M);
    y_P(1) = 0;
    y_P(end) = H;

end