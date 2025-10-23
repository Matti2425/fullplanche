function [vertices, faces] = generate_box(length, width, height)
% generate_box generates the vertices and faces for a rectangular box.
%
% Syntax:
%   [vertices, faces] = generate_box(length, width, height)
%
% Inputs:
%   length - Length along the x-axis (m)
%   width  - Width along the y-axis (m)
%   height - Height along the z-axis (m)
%
% Outputs:
%   vertices - Nx3 array of vertex coordinates
%   faces    - Mx4 array of face indices

% Define the vertices of the box centered at the origin
lx = length / 2;
ly = width / 2;
lz = height / 2;
vertices = [-lx, -ly, -lz;
             lx, -ly, -lz;
             lx,  ly, -lz;
            -lx,  ly, -lz;
            -lx, -ly,  lz;
             lx, -ly,  lz;
             lx,  ly,  lz;
            -lx,  ly,  lz];

% Define the faces of the box
faces = [1 2 3 4;  % Bottom face
         5 6 7 8;  % Top face
         1 2 6 5;  % Side face
         2 3 7 6;  % Side face
         3 4 8 7;  % Side face
         4 1 5 8]; % Side face

end
