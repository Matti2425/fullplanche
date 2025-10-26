function [psi, iter_outer, iter_inner] = solveSystem(a_n, a_s, a_e, a_w, a_p, a_p_base, b_p, ...
                                                      psi_in, psi_prev, Numerical)
% SOLVESYSTEM Solves the convection-diffusion equation using nested iteration
%
% This function implements a two-level iterative solver:
% - Outer loop: Time-stepping for transient behavior
% - Inner loop: Gauss-Seidel iteration for spatial discretization
%
% INPUTS:
%   a_n, a_s, a_e, a_w - Neighbor coefficients
%   a_p                - Central coefficient
%   a_p_base           - Transient term coefficient
%   b_p                - Source term (boundary conditions)
%   phi_in             - Initial guess for phi
%   phi_prev           - Previous time step solution
%   Numerical          - Structure with N, M, tolerance, maxIt, time_step
%
% OUTPUTS:
%   phi        - Converged solution
%   iter_outer - Number of outer (time) iterations
%   iter_inner - Total number of inner (Gauss-Seidel) iterations
%--------------------------------------------------------------------------

    N = Numerical.N;
    M = Numerical.M;
    tolerance = Numerical.tolerance;
    maxIt = Numerical.maxIt;
    time_step = Numerical.time_step;
    
    psi = psi_in;
    iter_inner = 0;
    
    % Outer loop: Time-stepping
    for iter_outer = 1:maxIt
        % Update source term with transient contribution
        b_p(2:N+1, 2:M+1) = a_p_base(2:N+1, 2:M+1) .* psi_prev(2:N+1, 2:M+1);
        
        % Inner loop: Gauss-Seidel iteration
        for m = 1:maxIt
            psi_old = psi;
            
            % Gauss-Seidel update (vectorized for interior nodes)
            for j = 2:M+1
                for i = 2:N+1
                    psi(i, j) = (a_e(i, j) * psi(i+1, j) + ...
                                 a_w(i, j) * psi(i-1, j) + ...
                                 a_n(i, j) * psi(i, j+1) + ...
                                 a_s(i, j) * psi(i, j-1) + ...
                                 b_p(i, j)) / a_p(i, j);
                end
            end
            
            % Boundary updates (vectorized)
            % Top boundary (j = M+2)
            i = 2:N+1;
            psi(i, M+2) = (a_e(i, M+2) .* psi(i+1, M+2) + ...
                          a_w(i, M+2) .* psi(i-1, M+2) + ...
                          a_s(i, M+2) .* psi(i, M+1) + ...
                          b_p(i, M+2)) ./ a_p(i, M+2);
            
            % Bottom boundary (j = 1)
            psi(i, 1) = (a_e(i, 1) .* psi(i+1, 1) + ...
                        a_w(i, 1) .* psi(i-1, 1) + ...
                        a_n(i, 1) .* psi(i, 2) + ...
                        b_p(i, 1)) ./ a_p(i, 1);
            
            % Left boundary (i = 1)
            j = 2:M+1;
            psi(1, j) = (a_e(1, j) .* psi(2, j) + ...
                        a_n(1, j) .* psi(1, j+1) + ...
                        a_s(1, j) .* psi(1, j-1) + ...
                        b_p(1, j)) ./ a_p(1, j);
            
            % Right boundary (i = N+2)
            psi(N+2, j) = (a_w(N+2, j) .* psi(N+1, j) + ...
                          a_n(N+2, j) .* psi(N+2, j+1) + ...
                          a_s(N+2, j) .* psi(N+2, j-1) + ...
                          b_p(N+2, j)) ./ a_p(N+2, j);
            
            % Check inner convergence
            err_phi = max(abs(psi - psi_old), [], 'all');
            iter_inner = iter_inner + 1;
            
            if err_phi < tolerance
                break;
            end
        end
        
        % Check outer convergence (transient)
        err = abs(psi_prev - psi);
        error_value = max(err / time_step, [], 'all');
        psi_prev = psi;
        
        if error_value <= tolerance
            break;
        end
    end
    
    if iter_outer == maxIt
        warning('Maximum outer iterations reached without convergence');
    end

end