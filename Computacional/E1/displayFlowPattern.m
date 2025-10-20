function displayFlowPattern(streamFunc, gridX, gridY)
% DISPLAYFLOWPATTERN Plots streamlines using stream function contours
%
% Generates a streamline visualization by contouring the stream function
% values over a computational grid. Automatically adjusts contour levels
% based on the stream function range.
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
%   Adaptively creates streamline contours in black color
%--------------------------------------------------------------------------

    [meshX, meshY] = meshgrid(gridX, gridY);
    
    % Get the range of stream function values (excluding boundary padding)
    psi_min = min(streamFunc(:));
    psi_max = max(streamFunc(:));
    psi_range = psi_max - psi_min;
    
    % Create adaptive contour levels to ensure visibility
    if psi_range > 0
        % Generate 30 levels evenly spaced across the actual range
        levels = linspace(psi_min + 0.01*psi_range, psi_max - 0.01*psi_range, 30);
        contour(meshX, meshY, streamFunc, levels, 'k', 'LineWidth', 0.01);
    else
        % If range is zero or very small, create synthetic levels
        contour(meshX, meshY, streamFunc, 30, 'k', 'LineWidth', 0.01);
    end

end