function [coeff_east, coeff_west, coeff_north, coeff_south, coeff_center, rhs] = ...
         constructMatrix(streamFunc, psi_boundary, solver, domain, velocity, gridX, gridY, fluidMask, state)
% CONSTRUCTMATRIX Builds coefficient matrices for discretized stream equation
%
% Constructs the finite difference stencil coefficients for solving the
% Laplace equation governing potential flow using Gauss-Seidel iteration.
%
% INPUTS:
%   streamFunc  - Current stream function field
%   psi_boundary - Stream function value on obstacle surface
%   solver      - Structure with nx, ny
%   domain      - Structure with length and height
%   velocity    - Inlet velocity
%   gridX, gridY - Coordinate vectors
%   fluidMask   - Binary matrix for solid/fluid regions
%   state       - Initial state structure
%
% OUTPUTS:
%   coeff_east, coeff_west, coeff_north, coeff_south - Face coefficients
%   coeff_center - Central point coefficient
%   rhs - Right-hand side vector
%--------------------------------------------------------------------------

    nx = solver.nx;
    ny = solver.ny;
    domainLength = domain.length;
    domainHeight = domain.height;

    dx = domainLength / nx;
    dy = domainHeight / ny;

    % Allocate coefficient matrices
    coeff_east = zeros(ny+2, nx+2);
    coeff_west = zeros(ny+2, nx+2);
    coeff_north = zeros(ny+2, nx+2);
    coeff_south = zeros(ny+2, nx+2);
    coeff_center = zeros(ny+2, nx+2);
    rhs = zeros(ny+2, nx+2);

    psi_bottom = state.psi_ref;

    % Interior nodes - construct finite difference stencil
    for j = 2:ny+1
        for i = 2:nx+1
            coeff_east(j, i) = dy / (gridX(i+1) - gridX(i));
            coeff_west(j, i) = dy / (gridX(i) - gridX(i-1));
            coeff_north(j, i) = dx / (gridY(j+1) - gridY(j));
            coeff_south(j, i) = dx / (gridY(j) - gridY(j-1));
            coeff_center(j, i) = coeff_east(j, i) + coeff_west(j, i) + ...
                                coeff_north(j, i) + coeff_south(j, i);
            rhs(j, i) = 0;
        end
    end

    % Bottom boundary (j=1) - Dirichlet condition
    for i = 1:nx+2
        coeff_east(1, i) = 0;
        coeff_west(1, i) = 0;
        coeff_north(1, i) = 0;
        coeff_south(1, i) = 0;
        coeff_center(1, i) = 1;
        rhs(1, i) = psi_bottom;
    end

    % Top boundary (j=ny+2) - Dirichlet condition
    for i = 1:nx+2
        coeff_east(ny+2, i) = 0;
        coeff_west(ny+2, i) = 0;
        coeff_north(ny+2, i) = 0;
        coeff_south(ny+2, i) = 0;
        coeff_center(ny+2, i) = 1;
        rhs(ny+2, i) = psi_bottom + velocity * domainHeight;
    end

    % Inlet boundary (i=1) - Dirichlet condition
    for j = 2:ny+1
        coeff_east(j, 1) = 0;
        coeff_west(j, 1) = 0;
        coeff_north(j, 1) = 0;
        coeff_south(j, 1) = 0;
        coeff_center(j, 1) = 1;
        rhs(j, 1) = psi_bottom + velocity * gridY(j);
    end

    % Outlet boundary (i=nx+2) - Neumann-like condition
    for j = 2:ny+1
        coeff_east(j, nx+2) = 0;
        coeff_west(j, nx+2) = 1;
        coeff_north(j, nx+2) = 0;
        coeff_south(j, nx+2) = 0;
        coeff_center(j, nx+2) = 1;
        rhs(j, nx+2) = 0;
    end

    % Obstacle surface - Dirichlet condition
    for j = 2:ny+1
        for i = 2:nx+1
            if fluidMask(j, i) == 0
                coeff_east(j, i) = 0;
                coeff_west(j, i) = 0;
                coeff_north(j, i) = 0;
                coeff_south(j, i) = 0;
                coeff_center(j, i) = 1;
                rhs(j, i) = psi_boundary;
            end
        end
    end

end