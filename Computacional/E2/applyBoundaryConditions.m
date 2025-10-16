function [a_n, a_s, a_e, a_w, a_p, b_p] = applyBoundaryConditions(a_n, a_s, a_e, a_w, a_p, a_p_base, b_p, ...
                                                                  surf_x, Numerical, diag, BoundaryConditions)
% APPLYBOUNDARYCONDITIONS Applies boundary conditions to coefficient matrices
%
% This function modifies the coefficient matrices to enforce boundary conditions
% for either the Smith-Hutton test case or the diagonal flow case.
%
% INPUTS:
%   a_n, a_s, a_e, a_w - Neighbor coefficients
%   a_p                - Central coefficient
%   a_p_base           - Transient term coefficient
%   b_p                - Source term
%   surf_x             - X-coordinates of mesh grid
%   Numerical          - Structure with N, M
%   diag               - Flag: 1 for diagonal, 0 for Smith-Hutton
%   BoundaryConditions - Structure with boundary values
%
% OUTPUTS:
%   a_n, a_s, a_e, a_w - Modified neighbor coefficients
%   a_p                - Modified central coefficient
%   b_p                - Modified source term
%--------------------------------------------------------------------------

    N = Numerical.N;
    M = Numerical.M;
    
    if diag == 1
        % Diagonal flow case: simple Dirichlet boundaries
        b_p(1, 1:M+2) = BoundaryConditions.phi_left;      % Left wall
        b_p(N+2, 1:M+2) = BoundaryConditions.phi_right;   % Right wall
        b_p(1:N+2, 1) = BoundaryConditions.phi_bottom;    % Bottom wall
        b_p(1:N+2, M+2) = BoundaryConditions.phi_top;     % Top wall
        
    else
        % Smith-Hutton test case boundaries
        % Left, right, and top boundaries
        b_p(1, 2:M+1) = 1 - tanh(10);
        b_p(N+2, 2:M+1) = 1 - tanh(10);
        b_p(2:N+1, M+2) = 1 - tanh(10);
        
        % Bottom boundary (varies with x position)
        for i = 2:N+1
            if surf_x(i, M+2) < 0
                % Left side of bottom boundary
                b_p(i, 1) = 1 + tanh(10 * (2 * surf_x(i, M+2) + 1));
            else
                % Right side of bottom boundary (Neumann condition)
                b_p(i, 1) = a_n(i, 1) + a_s(i, 1) + a_e(i, 1) + a_w(i, 1) + a_p_base(i, 1);
                a_n(i, 1) = 1;
            end
        end
    end

end