function [x_cyl, y_cyl, solid, theta] = cylinder_generator(x_center, y_center, R, M, N, x_P, y_P, solid)

% This function creates a set of points defining the perimeter of a cylinder 
% centered at (x_center, y_center) with a radius R. Additionally, it updates 
% the solid matrix to mark which grid points lie inside the cylinder.
%
% Inputs:
%   x_center, y_center - Coordinates of the cylinder center.
%   R - Radius of the cylinder.
%   M, N - Number of discretization points in the y and x directions.
%   x_P, y_P - Node points along the x and y directions.
%   solid - Predefined matrix indicating solid (0) or fluid (1) regions.
%
% Outputs:
%   x_cyl, y_cyl - Arrays containing the coordinates of the cylinder's perimeter.
%   solid - Updated matrix marking solid nodes inside the cylinder.
%--------------------------------------------------------------------------


theta = linspace(0, 2*pi, 10000); 
x_cyl = x_center + R * cos(theta);
y_cyl = y_center + R * sin(theta);


for j = 1:M+2
    for i = 1:N+2
        if (x_P(i) - x_center)^2 + (y_P(j) - y_center)^2 <= R^2
            solid(j, i) = 0;  
        end
    end
end
end


