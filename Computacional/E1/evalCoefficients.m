function [a_E, a_W, a_N, a_S, b_P, a_P] = evalCoefficients(phi, phi_body, NumericalParameters,Geometry, Vin, X, Y, S)

% This function constructs the coefficient matrices required for solving the discretized 
% form of the equation governing the stream function in a potential flow.
%
% INPUTS:
%   psi_start - Initial value of the stream function at the bot boundary.
%   psi_body  - Stream function value imposed on solid boundaries.
%   M         - Number of grid points in the vertical direction.
%   N         - Number of grid points in the horizontal direction.
%   L         - Length of the computational domain.
%   H         - Height of the computational domain.
%   v_in      - Free-stream velocity.
%   y_P, x_P  - Mesh grid points in the y and x directions.
%   solid     - Matrix indicating solid regions (0 for solid, 1 for fluid).
%
% OUTPUTS:
%   a_E, a_W, a_N, a_S - Discretization coefficients for east, west, north, and south neighbors.
%   a_P                 - Central coefficient of the discretized equation.
%   b_P                 - Right-hand side term of the equation.
%--------------------------------------------------------------------------

% Local variables
N = NumericalParameters.N;
M = NumericalParameters.M;

% Geometry increments
dx = (Geometry.L/N);   
dy = (Geometry.H/M); 

% Allocating space
a_E = zeros(N+2, M+2); 
a_W = zeros(N+2, M+2);
a_N = zeros(N+2, M+2);
a_S = zeros(N+2, M+2);
a_P = zeros(N+2, M+2);
b_P = zeros(N+2, M+2);

% Pre-defining values
phi_bot=phi;

for i = 2:N+1
    for j = 2:M+1
        a_E(i, j) = 1 * dy / (X(i+1)-X(i));
        a_W(i, j) = 1 * dy / (X(i)-X(i-1));
        a_N(i, j) = 1 * dx / (Y(j+1)-Y(j));
        a_S(i, j) = 1 * dx / (Y(j)-Y(j-1));
        a_P(i, j) = a_E(i, j) + a_W(i, j) + a_N(i, j) + a_S(i, j);
        b_P(i, j) = 0;
    end
end

for i = 1:N+2
    % Bot row
    a_E(1, i) = 0;
    a_W(1, i) = 0;
    a_N(1, i) = 0;
    a_S(1, i) = 0;
    a_P(1, i) = 1;
    b_P(1, i) = phi_bot;

    % Top row
    a_E(M+2, i) = 0;
    a_W(M+2, i) = 0;
    a_N(M+2, i) = 0;
    a_S(M+2, i) = 0;
    a_P(M+2, i) = 1;
    b_P(M+2, i) = phi_bot + v_in * H;
end

for j = 2:M+1
    % Inlet
    a_E(j, 1) = 0;
    a_W(j, 1) = 0;
    a_N(j, 1) = 0;
    a_S(j, 1) = 0;
    a_P(j, 1) = 1;
    b_P(j, 1) = phi_bot + v_in * Y(j);
    
    % Outlet
    a_E(j, N+2) = 0;
    a_W(j, N+2) = 1;
    a_N(j, N+2) = 0;
    a_S(j, N+2) = 0;
    a_P(j, N+2) = 1;
    b_P(j, N+2) = 0;

end

% Solid nodes
a_E(S) = 0;
a_W(S) = 0;
a_N(S) = 0;
a_S(S) = 0;
a_P(S) = 1;
b_P(S) = phi_body;

end

