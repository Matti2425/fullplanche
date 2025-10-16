function [vertices, faces] = generate_cube(side)
% generate_cube generates the vertices and faces for a cube.
%
% Syntax:
%   [vertices, faces] = generate_cube(side)
%
% Input:
%   side - Side length of the cube (m)
%
% Outputs:
%   vertices - Nx3 array of vertex coordinates
%   faces    - Mx4 array of face indices

% Define the vertices of the cube centered at the origin
vertices = [-1, -1, -1;
             1, -1, -1;
             1,  1, -1;
            -1,  1, -1;
            -1, -1,  1;
             1, -1,  1;
             1,  1,  1;
            -1,  1,  1];

% Scale the cube
halfSide = side / 2;
vertices = halfSide * vertices;

% Define the faces of the cube
faces = [1 2 3 4;  % Bottom face
         5 6 7 8;  % Top face
         1 2 6 5;  % Side face
         2 3 7 6;  % Side face
         3 4 8 7;  % Side face
         4 1 5 8]; % Side face

end
