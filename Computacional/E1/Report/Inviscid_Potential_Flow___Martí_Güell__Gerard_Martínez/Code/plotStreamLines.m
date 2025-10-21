function plotStreamLines(psi, x_P, y_P)
% PLOTSTREAMLINES Plots the streamlines of a flow field using contour lines.
%
% This function generates a streamline plot by contouring the stream function 
% `psi` over a structured grid defined by `x_P` and `y_P`. The function 
% represents streamlines that indicate the direction of the flow.
%
% Inputs:
%   psi  - Matrix of stream function values.
%   x_P  - Vector of x-coordinates for the grid points.
%   y_P  - Vector of y-coordinates for the grid points.
%
% The function uses `contour` to draw 30 evenly spaced streamlines with 
% black color and a thin line width. The contour lines represent the flow 
% paths in a 2D domain.
%
% Outputs:
%   A graphical representation of streamlines.
%--------------------------------------------------------------------------

[X, Y] = meshgrid(x_P, y_P);

contour(X, Y, psi, 30, 'k', 'LineWidth', 0.01);

end