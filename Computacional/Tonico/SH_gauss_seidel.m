function [Phi, iter, error] = SH_gauss_seidel(a_E, a_W, a_N, a_S, a_P, b, delta, M, N, alpha, x)

% Solves the 2D convection-diffusion equation using the
% Gauss-Seidel iterative method for a structured mesh.
%
% The algorithm updates the scalar field Phi at the internal nodes using
% the discrete equation coefficients and applies Dirichlet and Neumann
% boundary conditions at the inlet and outlet, respectively.
%
% Inputs:
%   - a_E, a_W, a_N, a_S : Discrete coefficients for each face
%   - a_P                : Central coefficient for each node
%   - b                  : Right-hand side (source term)
%   - delta              : Convergence tolerance
%   - M, N               : Number of cells in y and x directions
%   - alpha              : Parameter controlling the inlet boundary profile
%   - x                  : x-coordinate of cell centers
%
% Outputs:
%   - Phi   : Scalar field solution (2D matrix)
%   - iter  : Number of iterations performed
%   - error : Final error (residual) of the iterative method
%---------------------------------------------------------------------------

Phi=ones(M+2, N+2)*(1 - tanh(alpha));

max_iter = 100000;
error = Inf;

for iter = 1:max_iter
    phi_old = Phi;  % Guardamos el valor anterior para evaluar el error


        %Inlet
        for i=1:ceil((N+2)/2)
        Phi(1,i)= 1 + tanh (alpha*(2*x(i) + 1));
        end

    for j = 2:M+1
        for i = 2:N+1
            Phi(j, i) = (a_E(j, i) * Phi(j, i+1) + ...
                         a_W(j, i) * Phi(j, i-1) + ...
                         a_N(j, i) * Phi(j+1, i) + ...
                         a_S(j, i) * Phi(j-1, i) + ...
                         b(j, i)) / a_P(j, i);
        end
    end

        %Outlet
        for i=ceil((N+2)/2+1):N+2
        Phi(1,i)= Phi(2,i);
        end

    % Convergence criteria
    error = max(max(abs(Phi - phi_old)));

    if error < delta
        break;
    end
end

end
