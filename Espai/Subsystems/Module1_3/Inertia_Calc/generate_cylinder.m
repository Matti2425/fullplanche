function [X, Y, Z] = generate_cylinder(radius, height)
% generate_cylinder generates the mesh for a Z-aligned cylinder.
%
% Syntax:
%   [X, Y, Z] = generate_cylinder(radius, height)
%
% Inputs:
%   radius - Radius of the cylinder's base (m)
%   height - Height of the cylinder (m)
%
% Outputs:
%   X, Y, Z - Meshgrid coordinates of the cylinder surface

% Generate a cylinder along the z-axis
[Xc, Yc, Zc] = cylinder(radius, 20); % 20 divisions for smoothness
Zc = Zc * height - height / 2; % Center the cylinder at z=0

% Adjust the cylinder to the Z axis
X = Xc;
Y = Yc;
Z = Zc;

end
