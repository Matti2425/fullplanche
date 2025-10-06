function [psi, iter, error] = gauss_seidel(psi, psi_body, solid, a_E, a_W, a_N, a_S, a_P, b_P, v_in, y_P, M, N, H, delta)

% GAUSS_SEIDEL Solves the stream function equation using the Gauss-Seidel method.
%
% This function iteratively solves the discretized stream function equation 
% using the Gauss-Seidel method with Dirichlet boundary conditions and 
% internal conditions imposed by solid obstacles in the domain.
%
% Inputs:
%   psi       - Initial guess for the stream function.
%   psi_body  - Stream function value inside solid obstacles.
%   solid     - Binary matrix indicating solid (0) or fluid (1) regions.
%   a_E, a_W, a_N, a_S, a_P - Coefficients of the discretized equation.
%   b_P       - Right-hand side of the equation.
%   v_in      - Inlet velocity.
%   y_P       - Array of y-coordinates of the grid points.
%   M, N      - Number of grid points in the y and x directions.
%   H         - Channel height.
%   delta     - Convergence threshold for the iterative solver.
%
% Outputs:
%   psi       - Computed stream function field.
%   iter      - Number of iterations performed.
%   error     - Final error value.
%
% The solver iterates until convergence is achieved, applying boundary 
% conditions at the channel walls and enforcing the presence of solid 
% obstacles. The Gauss-Seidel method is used to update the stream function 
% values at each iteration. The algorithm stops when the maximum difference 
% between successive iterations is below the specified tolerance (delta).
%--------------------------------------------------------------------------

max_iter = 1000000;
error = Inf; 
for iter = 1:max_iter
    psi_old = psi; 
    % Channel enter (x = 0) 
    for j = 1:M+2
        psi(j, 1) = v_in * y_P(j);
    end
    % Channel exit (x = L) 
    for j = 1:M+2
        psi(j, N+2) = psi(j, N+1);
    end  
    % Inferior wall
    for i = 1:N+2
        psi(1, i) = 0;
    end
    % Superior wall
    for i = 1:N+2
        psi(M+2, i) = v_in * H;
    end
    % Solid
    for j = 1:M+2
        for i = 1:N+2
            if solid(j, i) == 0
                psi(j, i) = psi_body;  
            end
        end
    end

    for j = 2:M+1
        for i = 2:N+1
                psi(j, i) = (a_E(j, i) * psi(j, i+1) + a_W(j, i) * psi(j, i-1) + ...
                             a_N(j, i) * psi(j+1, i) + a_S(j, i) * psi(j-1, i) + b_P(j, i)) / a_P(j, i);
        end
    end
    error = max(max(abs(psi - psi_old)));
    if error < delta
        break;
    end
end
end

