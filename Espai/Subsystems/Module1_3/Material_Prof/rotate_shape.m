function [X_rot, Y_rot, Z_rot] = rotate_shape(X, Y, Z, R)
% rotate_shape applies a rotation matrix to the coordinates X, Y, Z.
%
% Syntax:
%   [X_rot, Y_rot, Z_rot] = rotate_shape(X, Y, Z, R)
%
% Inputs:
%   X, Y, Z - Coordinate matrices of the shape
%   R       - 3x3 rotation matrix
%
% Outputs:
%   X_rot, Y_rot, Z_rot - Rotated coordinate matrices

% Flatten the coordinate matrices
X_flat = X(:);
Y_flat = Y(:);
Z_flat = Z(:);

% Combine into a matrix
coords = [X_flat'; Y_flat'; Z_flat'];

% Apply rotation
coords_rot = R * coords;

% Reshape back to original shape
X_rot = reshape(coords_rot(1,:), size(X));
Y_rot = reshape(coords_rot(2,:), size(Y));
Z_rot = reshape(coords_rot(3,:), size(Z));

end
