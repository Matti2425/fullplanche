function [a_E, a_W, a_N, a_S, b_P, a_P] = evalCoefficients(psi, psi_body, NumericalParameters, Geometry, v_in, x_P, y_P, solid, Init)
% EVALCOEFFICIENTS Constructs coefficient matrices for discretized stream function equation
%
% This function constructs the coefficient matrices required for solving the discretized 
% form of the equation governing the stream function in a potential flow.
%
% INPUTS:
%   psi                 - Current stream function field
%   psi_body            - Stream function value imposed on solid boundaries
%   NumericalParameters - Structure with N, M
%   Geometry            - Structure with L, H
%   v_in                - Free-stream velocity
%   x_P, y_P            - 1D arrays of mesh grid points
%   solid               - Matrix indicating solid regions (0 for solid, 1 for fluid)
%   Init                - Structure with phi_start
%
% OUTPUTS:
%   a_E, a_W, a_N, a_S - Discretization coefficients for east, west, north, and south neighbors
%   a_P                 - Central coefficient of the discretized equation
%   b_P                 - Right-hand side term of the equation
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    L = Geometry.L;
    H = Geometry.H;

    dx = L / N;   
    dy = H / M; 

    % Allocate coefficient matrices (M+2, N+2)
    a_E = zeros(M+2, N+2); 
    a_W = zeros(M+2, N+2);
    a_N = zeros(M+2, N+2);
    a_S = zeros(M+2, N+2);
    a_P = zeros(M+2, N+2);
    b_P = zeros(M+2, N+2);

    psi_bottom = Init.phi_start;

    % Interior nodes - vectorized computation
    for j = 2:M+1
        for i = 2:N+1
            a_E(j, i) = dy / (x_P(i+1) - x_P(i));
            a_W(j, i) = dy / (x_P(i) - x_P(i-1));
            a_N(j, i) = dx / (y_P(j+1) - y_P(j));
            a_S(j, i) = dx / (y_P(j) - y_P(j-1));
            a_P(j, i) = a_E(j, i) + a_W(j, i) + a_N(j, i) + a_S(j, i);
            b_P(j, i) = 0;
        end
    end

    % Bottom boundary (j=1)
    for i = 1:N+2
        a_E(1, i) = 0;
        a_W(1, i) = 0;
        a_N(1, i) = 0;
        a_S(1, i) = 0;
        a_P(1, i) = 1;
        b_P(1, i) = psi_bottom;
    end

    % Top boundary (j=M+2)
    for i = 1:N+2
        a_E(M+2, i) = 0;
        a_W(M+2, i) = 0;
        a_N(M+2, i) = 0;
        a_S(M+2, i) = 0;
        a_P(M+2, i) = 1;
        b_P(M+2, i) = psi_bottom + v_in * H;
    end

    % Inlet boundary (i=1)
    for j = 2:M+1
        a_E(j, 1) = 0;
        a_W(j, 1) = 0;
        a_N(j, 1) = 0;
        a_S(j, 1) = 0;
        a_P(j, 1) = 1;
        b_P(j, 1) = psi_bottom + v_in * y_P(j);
    end

    % Outlet boundary (i=N+2)
    for j = 2:M+1
        a_E(j, N+2) = 0;
        a_W(j, N+2) = 1;
        a_N(j, N+2) = 0;
        a_S(j, N+2) = 0;
        a_P(j, N+2) = 1;
        b_P(j, N+2) = 0;
    end

    % Solid nodes
    for j = 2:M+1
        for i = 2:N+1
            if solid(j, i) == 0
                a_E(j, i) = 0;
                a_W(j, i) = 0;
                a_N(j, i) = 0;
                a_S(j, i) = 0;
                a_P(j, i) = 1;
                b_P(j, i) = psi_body;
            end
        end
    end

end

