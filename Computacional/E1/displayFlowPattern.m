function displayFlowPattern(streamFunc, gridX, gridY)
% DISPLAYFLOWPATTERN Plots streamlines using stream function contours
%
% Generates a streamline visualization by contouring the stream function
% values over a computational grid.
%
% Inputs:
%   streamFunc - Stream function values at grid points (ny+2, nx+2)
%   gridX      - X-coordinate vector
%   gridY      - Y-coordinate vector
%
% Outputs:
%   Graphical streamline contour plot
%
% Notes:
%   Creates 30 streamline contours in black color
%--------------------------------------------------------------------------

    [meshX, meshY] = meshgrid(gridX, gridY);
    contour(meshX, meshY, streamFunc, 30, 'k', 'LineWidth', 0.01);

end