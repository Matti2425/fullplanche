function [X, Y, Z] = generate_sphere(radius)
% generate_sphere generates the mesh for a sphere.
%
% Syntax:
%   [X, Y, Z] = generate_sphere(radius)
%
% Input:
%   radius - Radius of the sphere (m)
%
% Outputs:
%   X, Y, Z - Meshgrid coordinates of the sphere surface

% Generate sphere mesh with 20 divisions for smoothness
[X, Y, Z] = sphere(20);
X = X * radius;
Y = Y * radius;
Z = Z * radius;

end
