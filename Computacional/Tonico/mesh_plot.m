function mesh_plot(x_P, y_P, L, H, N, M)

% INPUTS:
%   x_P - X-coordinates of the mesh nodes.
%   y_P - Y-coordinates of the mesh nodes.
%   L   - Length of the computational domain.
%   H   - Height of the computational domain.
%   N   - Number of nodes in the X-direction.
%   M   - Number of nodes in the Y-direction.
%
% OUTPUT:
%   A 2D plot displaying the mesh structure, including node positions, 
%   cell divisions, and domain boundaries.
%--------------------------------------------------------------------------

dx = L/N;   
dy = H/M; 

[X, Y] = meshgrid(x_P, y_P);
    
% Plot the nodes
plot(X(:), Y(:), 'b.', 'MarkerSize', 0.5, 'MarkerFaceColor', 'b'); 
axis equal
axis tight
hold on;

for i = 2:N
    plot([x_P(i)+0.5*dx, x_P(i)+0.5*dx], [0, H], 'k--', 'LineWidth', 0.1);
end
    
% Horizontal lines
for j = 2:M
    plot([0, L], [y_P(j)+0.5*dy, y_P(j)+0.5*dy], 'k--', 'LineWidth', 0.1);
end

%Channel walls
plot([0, L], [0, 0], 'k-', 'LineWidth', 1); % Línea inferior en negro y más gruesa
plot([0, L], [H, H], 'k-', 'LineWidth', 1); % Línea superior en negro y más gruesa
plot([0, 0], [0, H], 'k--', 'LineWidth', 1); % Línea izquierda
plot([L, L], [0, H], 'k--', 'LineWidth', 1); % Línea derecha
 axis equal
 axis tight
title('Malla de cálculo con nodos en el centro de las celdas');
xlabel('x [m]'); ylabel('y [m]'); grid off;

end

