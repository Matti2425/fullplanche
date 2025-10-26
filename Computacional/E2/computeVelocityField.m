function [v_x, v_y] = computeVelocityField(surf_x, surf_y, Fluid, diag)
% COMPUTEVELOCITYFIELD Computes the velocity field for the convection-diffusion problem
%
% This function calculates the velocity components based on the problem type:
% - Smith-Hutton test case: v_x = 2*y*(1-x^2), v_y = -2*x*(1-y^2)
% - Diagonal flow: v_x = v_0*cos(45°), v_y = v_0*sin(45°)
%
% INPUTS:
%   surf_x, surf_y - Meshgrid coordinates
%   Fluid          - Structure containing fluid properties (v_0 for diagonal case)
%   diag           - Flag: 1 for diagonal flow, 0 for Smith-Hutton
%
% OUTPUTS:
%   v_x - Velocity component in x-direction
%   v_y - Velocity component in y-direction
%--------------------------------------------------------------------------

    if diag
        % Diagonal flow case: uniform velocity at 45 degrees
        v_x = Fluid.v_0 * cosd(45) * ones(size(surf_x));
        v_y = Fluid.v_0 * sind(45) * ones(size(surf_y));
    else
        % Smith-Hutton test case: analytical velocity field
        v_x = 2 * surf_y .* (1 - surf_x.^2);
        v_y = -2 * surf_x .* (1 - surf_y.^2);
    end

end