function [a_E, a_W, a_N, a_S, a_P, rhs] = ...
         constructMatrix(psi, psi_boundary, solver, domain, velocity, gridX, gridY, fluidMask, state)
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
    a_E = zeros(ny+2, nx+2);
    a_W = zeros(ny+2, nx+2);
    a_N = zeros(ny+2, nx+2);
    a_S = zeros(ny+2, nx+2);
    a_P = zeros(ny+2, nx+2);
    rhs = zeros(ny+2, nx+2);

    psi_bottom = state.psi_ref;

    % Interior nodes - construct finite difference stencil
    for j = 2:ny+1
        for i = 2:nx+1
            a_E(j, i) = dy / (gridX(i+1) - gridX(i));
            a_W(j, i) = dy / (gridX(i) - gridX(i-1));
            a_N(j, i) = dx / (gridY(j+1) - gridY(j));
            a_S(j, i) = dx / (gridY(j) - gridY(j-1));
            a_P(j, i) = a_E(j, i) + a_W(j, i) + ...
                                a_N(j, i) + a_S(j, i);
            rhs(j, i) = 0;
        end
    end

    % Bottom boundary (j=1) - Dirichlet condition
    for i = 1:nx+2
        a_E(1, i) = 0;
        a_W(1, i) = 0;
        a_N(1, i) = 0;
        a_S(1, i) = 0;
        a_P(1, i) = 1;
        rhs(1, i) = psi_bottom;
    end

    % Top boundary (j=ny+2) - Dirichlet condition
    for i = 1:nx+2
        a_E(ny+2, i) = 0;
        a_W(ny+2, i) = 0;
        a_N(ny+2, i) = 0;
        a_S(ny+2, i) = 0;
        a_P(ny+2, i) = 1;
        rhs(ny+2, i) = psi_bottom + velocity * domainHeight;
    end

    % Inlet boundary (i=1) - Dirichlet condition
    for j = 2:ny+1
        a_E(j, 1) = 0;
        a_W(j, 1) = 0;
        a_N(j, 1) = 0;
        a_S(j, 1) = 0;
        a_P(j, 1) = 1;
        rhs(j, 1) = psi_bottom + velocity * gridY(j);
    end

    % Outlet boundary (i=nx+2) - Neumann-like condition
    for j = 2:ny+1
        a_E(j, nx+2) = 0;
        a_W(j, nx+2) = 1;
        a_N(j, nx+2) = 0;
        a_S(j, nx+2) = 0;
        a_P(j, nx+2) = 1;
        rhs(j, nx+2) = 0;
    end

    % Obstacle surface - Dirichlet condition
    for j = 2:ny+1
        for i = 2:nx+1
            if fluidMask(j, i) == 0
                a_E(j, i) = 0;
                a_W(j, i) = 0;
                a_N(j, i) = 0;
                a_S(j, i) = 0;
                a_P(j, i) = 1;
                rhs(j, i) = psi_boundary;
            end
        end
    end

end