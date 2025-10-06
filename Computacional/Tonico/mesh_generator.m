function [x_P,y_P] = mesh_generator(N, M, L, H)

% MESH_GENERATOR Creates a structured Cartesian mesh for a 2D domain.
%
% This function generates a uniform computational grid for a domain 
% of dimensions L × H, discretized with N points in the x-direction 
% and M points in the y-direction.
%
% Inputs:
%   N, M - Number of grid points in the x and y directions, respectively.
%   L, H - Length and height of the domain.
%
% Outputs:
%   x_P - Array containing the x-coordinates of the cell centers and boundaries.
%   y_P - Array containing the y-coordinates of the cell centers and boundaries.
%
% The function initializes the mesh by defining the x and y coordinates
% of the grid points. It places the first and last points at the domain 
% boundaries, ensuring proper discretization.

x_P=zeros(1, N+2);
x_P(2:end-1) = (0.5:1:N-0.5) * (L/N);
x_P(1) = 0;
x_P(end) = L;

y_P=zeros(1, M+2);
y_P(2:end-1) = (0.5:1:M-0.5) * (H/M);
y_P(1) = 0;
y_P(end) = H;

end

