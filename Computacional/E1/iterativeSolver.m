function [streamFunc, numIter, residual] = iterativeSolver(streamFunc, psi_boundary, fluidMask, ...
                                                           coeff_east, coeff_west, coeff_north, ...
                                                           coeff_south, coeff_center, rhs, ...
                                                           velocity, gridY, solver, domain)
% ITERATIVESOLVER Solves Laplace equation using Gauss-Seidel iteration
%
% Implements the Gauss-Seidel iterative method with Dirichlet boundary
% conditions and internal solid region constraints.
%
% Inputs:
%   streamFunc  - Initial stream function field (guess)
%   psi_boundary - Stream function on obstacle surface
%   fluidMask   - Binary matrix for solid/fluid regions
%   coeff_*     - Finite difference coefficients
%   rhs         - Right-hand side vector
%   velocity    - Inlet velocity
%   gridY       - Y-coordinate vector
%   solver      - Structure with tolerance, max_iterations, nx, ny
%   domain      - Structure with height
%
% Outputs:
%   streamFunc - Converged stream function field
%   numIter    - Number of iterations performed
%   residual   - Final convergence error
%--------------------------------------------------------------------------

    nx = solver.nx;
    ny = solver.ny;
    domainHeight = domain.height;
    tolerance = solver.tolerance;
    maxIterations = solver.max_iterations;

    residual = Inf;

    for numIter = 1:maxIterations
        streamFuncOld = streamFunc;

        % Apply inlet boundary condition
        for j = 1:ny+2
            streamFunc(j, 1) = velocity * gridY(j);
        end

        % Apply outlet boundary condition
        for j = 1:ny+2
            streamFunc(j, nx+2) = streamFunc(j, nx+1);
        end

        % Apply bottom wall boundary condition
        for i = 1:nx+2
            streamFunc(1, i) = 0;
        end

        % Apply top wall boundary condition
        for i = 1:nx+2
            streamFunc(ny+2, i) = velocity * domainHeight;
        end

        % Apply obstacle boundary condition
        for j = 1:ny+2
            for i = 1:nx+2
                if fluidMask(j, i) == 0
                    streamFunc(j, i) = psi_boundary;
                end
            end
        end

        % Gauss-Seidel iteration over interior nodes
        for j = 2:ny+1
            for i = 2:nx+1
                streamFunc(j, i) = (coeff_east(j, i) * streamFunc(j, i+1) + ...
                                   coeff_west(j, i) * streamFunc(j, i-1) + ...
                                   coeff_north(j, i) * streamFunc(j+1, i) + ...
                                   coeff_south(j, i) * streamFunc(j-1, i) + ...
                                   rhs(j, i)) / coeff_center(j, i);
            end
        end

        % Check convergence criterion
        residual = max(max(abs(streamFunc - streamFuncOld)));
        if residual < tolerance
            break;
        end
    end

end