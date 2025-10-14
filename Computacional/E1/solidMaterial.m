function [x_cyl, y_cyl, solid, theta] = solidMaterial(Geometry, x_P, y_P, NumericalParameters)
% CYLINDER_GENERATOR Creates a cylinder in the computational domain
%
% This function creates a set of points defining the perimeter of a cylinder 
% centered at ObjectCenter with radius D/2. Additionally, it updates 
% the solid matrix to mark which grid points lie inside the cylinder.
%
% Inputs:
%   Geometry            - Structure containing ObjectCenter [x,y] and D (diameter)
%   x_P, y_P            - 1D arrays of grid coordinates
%   NumericalParameters - Structure containing N and M
%
% Outputs:
%   x_cyl, y_cyl - Arrays containing the coordinates of the cylinder's perimeter
%   solid        - Matrix marking solid (0) or fluid (1) regions
%   theta        - Angular positions around cylinder
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    x_center = Geometry.ObjectCenter(1);
    y_center = Geometry.ObjectCenter(2);
    R = Geometry.D / 2;  % Radius

    % Initialize solid matrix (1 = fluid, 0 = solid)
    solid = ones(M+2, N+2);

    % Generate cylinder perimeter
    theta = linspace(0, 2*pi, 10000);
    x_cyl = x_center + R * cos(theta);
    y_cyl = y_center + R * sin(theta);

    % Mark solid nodes inside the cylinder
    for j = 1:M+2
        for i = 1:N+2
            if (x_P(i) - x_center)^2 + (y_P(j) - y_center)^2 <= R^2
                solid(j, i) = 0;
            end
        end
    end

end