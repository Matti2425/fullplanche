function [a_E, a_W, a_N, a_S, a_P, b, AP_e, AP_w, AP_n, AP_s] = compute_coefficients(u, v, rho, Gamma, x, y, dx, dy, N, M, scheme)

% DESCRIPTION:
%   This function computes the coefficients of the discretized 
%   convection-diffusion. It supports multiple numerical schemes
%   (UDS, CDS, HDS, EDS, PLDS) to approximate the convective fluxes
%   based on the local Peclet number.
%
%   The function returns the coefficients of each cell face (a_E, a_W,
%   a_N, a_S), the central coefficient (a_P), the right-hand side term b,
%   and the values A(|P|) used in the numerical schemes.
%
% INPUTS:
%   u, v      - Velocity components in x and y directions (2D matrices)
%   rho       - Density field (assumed constant or given as 2D matrix)
%   Gamma     - Diffusivity field (2D matrix)
%   x, y      - Node coordinates in x and y directions
%   dx, dy    - Grid spacing in x and y directions
%   N, M      - Number of cells in x and y
%   scheme    - Integer indicating the numerical scheme:
%                1 - UDS (Upwind)
%                2 - CDS (Central)
%                3 - HDS (Hybrid)
%                4 - EDS (Exponential)
%                5 - PLDS (Power-law)
%
% OUTPUTS:
%   a_E, a_W, a_N, a_S, a_P, b - Discretization coefficients              
%   AP_e, AP_w, AP_n, AP_s - A(|P|) values for each face
%--------------------------------------------------------------------------

% Matrix preallocation 
F_e = zeros(M+2, N+2);
F_w = zeros(M+2, N+2);
F_n = zeros(M+2, N+2);
F_s = zeros(M+2, N+2);

D_e = zeros(M+2, N+2);
D_w = zeros(M+2, N+2);
D_n = zeros(M+2, N+2);
D_s = zeros(M+2, N+2);

a_E = zeros(M+2, N+2);
a_W = zeros(M+2, N+2);
a_N = zeros(M+2, N+2);
a_S = zeros(M+2, N+2);
a_P = zeros(M+2, N+2);
b   = zeros(M+2, N+2);

for i = 2:N+1
    for j = 2:M+1

        % Flujos convectivos
        F_e(j,i) = 0.5 * (rho(j,i+1)*u(j,i+1) + rho(j,i)*u(j,i)) * dy;
        F_w(j,i) = 0.5 * (rho(j,i)*u(j,i) + rho(j,i-1)*u(j,i-1)) * dy;
        F_n(j,i) = 0.5 * (rho(j+1,i)*v(j+1,i) + rho(j,i)*v(j,i)) * dx;
        F_s(j,i) = 0.5 * (rho(j,i)*v(j,i) + rho(j-1,i)*v(j-1,i)) * dx;

        % Coeficientes difusivos
        D_e(j,i) = (0.5 * (Gamma(j,i+1) + Gamma(j,i)) * dy) / (x(i+1) - x(i));
        D_w(j,i) = (0.5 * (Gamma(j,i) + Gamma(j,i-1)) * dy) / (x(i) - x(i-1));
        D_n(j,i) = (0.5 * (Gamma(j+1,i) + Gamma(j,i)) * dx) / (y(j+1) - y(j));
        D_s(j,i) = (0.5 * (Gamma(j,i) + Gamma(j-1,i)) * dx) / (y(j) - y(j-1));

        % Peclet number for each face
        P_e = F_e(j,i) / D_e(j,i);
        P_w = F_w(j,i) / D_w(j,i);
        P_n = F_n(j,i) / D_n(j,i);
        P_s = F_s(j,i) / D_s(j,i);

        % A(|P|) according to the numerical scheme
        switch scheme
            case 1  % UDS
                AP_e = 1;
                AP_w = 1;
                AP_n = 1;
                AP_s = 1;

            case 2  % CDS
                AP_e = 1 - 0.5 * abs(P_e);
                AP_w = 1 - 0.5 * abs(P_w);
                AP_n = 1 - 0.5 * abs(P_n);
                AP_s = 1 - 0.5 * abs(P_s);

            case 3  % HDS
                AP_e = max(0, 1 - 0.5 * abs(P_e));
                AP_w = max(0, 1 - 0.5 * abs(P_w));
                AP_n = max(0, 1 - 0.5 * abs(P_n));
                AP_s = max(0, 1 - 0.5 * abs(P_s));

            case 4  % EDS
                AP_e = abs(P_e) / (exp(abs(P_e)) - 1 + 1e-10);
                AP_w = abs(P_w) / (exp(abs(P_w)) - 1 + 1e-10);
                AP_n = abs(P_n) / (exp(abs(P_n)) - 1 + 1e-10);
                AP_s = abs(P_s) / (exp(abs(P_s)) - 1 + 1e-10);

            case 5  % PLDS
                AP_e = max(0, (1 - 0.5 * abs(P_e))^5);
                AP_w = max(0, (1 - 0.5 * abs(P_w))^5);
                AP_n = max(0, (1 - 0.5 * abs(P_n))^5);
                AP_s = max(0, (1 - 0.5 * abs(P_s))^5);

            otherwise
                error('Esquema numérico no reconocido');
        end

        % Compute a_i coefficients
        a_E(j,i) = D_e(j,i) * AP_e + max(-F_e(j,i), 0);
        a_W(j,i) = D_w(j,i) * AP_w + max( F_w(j,i), 0);
        a_N(j,i) = D_n(j,i) * AP_n + max(-F_n(j,i), 0);
        a_S(j,i) = D_s(j,i) * AP_s + max( F_s(j,i), 0);

        a_P(j,i) = a_E(j,i) + a_W(j,i) + a_N(j,i) + a_S(j,i);% Steady-state without source term
        b(j,i)   = 0; % Steady-state without source term

    end
end
end
