function [vertices, faces] = generate_tri_prism(length, width, height)
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
lx = length;
ly = width;
lz = height / 2;
vertices = [-lx*1/3, -ly*1/3, -lz;
             lx*2/3, -ly*1/3, -lz;
            -lx*1/3,  ly*2/3, -lz;
            -lx*1/3, -ly*1/3,  lz;
             lx*2/3, -ly*1/3,  lz;
            -lx*1/3,  ly*2/3,  lz];

% Define the faces of the box
faces = [1 2 3 1;  % Bottom face
         4 5 6 1;  % Top face
         1 2 4 5;  % Side face
         1 3 4 6;  % Side face
         2 3 5 6]; % Side face

end
