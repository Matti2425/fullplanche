function [psi, iter, error] = solveSystem(psi, psi_body, solid, a_E, a_W, a_N, a_S, a_P, b_P, v_in, y_P, NumericalParameters, Geometry)
% SOLVESYSTEM Solves the stream function equation using the Gauss-Seidel method
%
% This function iteratively solves the discretized stream function equation 
% using the Gauss-Seidel method with Dirichlet boundary conditions and 
% internal conditions imposed by solid obstacles in the domain.
%
% Inputs:
%   psi                 - Initial guess for the stream function (M+2, N+2)
%   psi_body            - Stream function value inside solid obstacles
%   solid               - Binary matrix indicating solid (0) or fluid (1) regions
%   a_E, a_W, a_N, a_S, a_P - Coefficients of the discretized equation
%   b_P                 - Right-hand side of the equation
%   v_in                - Inlet velocity
%   y_P                 - 1D array of y-coordinates of the grid points
%   NumericalParameters - Structure with N, M, i_max, d (tolerance)
%   Geometry            - Structure with H (height)
%
% Outputs:
%   psi   - Computed stream function field (M+2, N+2)
%   iter  - Number of iterations performed
%   error - Final error value
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    H = Geometry.H;
    delta = NumericalParameters.d;
    max_iter = NumericalParameters.i_max;

    error = Inf; 
    
    for iter = 1:max_iter
        psi_old = psi; 
        
        % Channel inlet (i = 1)
        for j = 1:M+2
            psi(j, 1) = v_in * y_P(j);
        end
        
        % Channel outlet (i = N+2)
        for j = 1:M+2
            psi(j, N+2) = psi(j, N+1);
        end
        
        % Bottom wall (j = 1)
        for i = 1:N+2
            psi(1, i) = 0;
        end
        
        % Top wall (j = M+2)
        for i = 1:N+2
            psi(M+2, i) = v_in * H;
        end
        
        % Solid nodes
        for j = 1:M+2
            for i = 1:N+2
                if solid(j, i) == 0
                    psi(j, i) = psi_body;
                end
            end
        end

        % Gauss-Seidel iteration for interior nodes
        for j = 2:M+1
            for i = 2:N+1
                psi(j, i) = (a_E(j, i) * psi(j, i+1) + a_W(j, i) * psi(j, i-1) + ...
                             a_N(j, i) * psi(j+1, i) + a_S(j, i) * psi(j-1, i) + b_P(j, i)) / a_P(j, i);
            end
        end
        
        % Check convergence
        error = max(max(abs(psi - psi_old)));
        if error < delta
            break;
        end
    end
    
end