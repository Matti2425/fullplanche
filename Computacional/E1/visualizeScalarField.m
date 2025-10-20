function visualizeScalarField(fieldLabel, fieldData, nx, ny, dx, dy)
% VISUALIZESCALARFIELD Creates a color-mapped visualization of a scalar field
%
% Generates a graphical representation of a 2D scalar field using
% interpolated shading and color mapping.
%
% Inputs:
%   fieldLabel - Title/label for the field (string)
%   fieldData  - 2D matrix of scalar values (ny+2, nx+2)
%   nx         - Number of grid points in x-direction
%   ny         - Number of grid points in y-direction
%   dx         - Grid spacing in x-direction
%   dy         - Grid spacing in y-direction
%
% No explicit outputs (creates figure)
%--------------------------------------------------------------------------

    % Display field with interpolated shading
    imagesc((0.5:nx-0.5) * dx, (0.5:ny-0.5) * dy, fieldData(2:end-1, 2:end-1));
    axis xy equal tight;
    cb = colorbar();
    cb.Label.String = fieldLabel;
    shading interp;
    colormap jet;

    % Configure axes and labels
    xlabel('x [m]');
    ylabel('y [m]');
    grid off;

    % Set title
    if isempty(fieldLabel)
        title('Scalar Field', 'FontWeight', 'bold', 'FontSize', 12);
    else
        title(['Field: ', fieldLabel], 'FontWeight', 'bold', 'FontSize', 12);
    end

end