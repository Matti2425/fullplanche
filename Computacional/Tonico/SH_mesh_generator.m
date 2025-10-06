function [x_P,y_P] = SH_mesh_generator(N, M, L, H)

% Generates the coordinates of cell centers (P nodes) in a 2D structured mesh
% over a domain of size L × H, with N divisions along x and M along y.
% Includes ghost cells at the boundaries.
%
% Inputs:
%   - N, M : Number of divisions in x and y directions
%   - L, H : Physical dimensions of the domain
%
% Outputs:
%   - x_P, y_P : Vectors with the coordinates of the cell centers
%--------------------------------------------------------------------------

x_P=zeros(1, N+2);
x_P(2:end-1) = (0.5-N/2:1:N/2-0.5) * (L/N);
x_P(1) = -L/2;
x_P(end) = L/2;

y_P=zeros(1, M+2);
y_P(2:end-1) = (0.5:1:M-0.5) * (H/M);
y_P(1) = 0;
y_P(end) = H;

end

