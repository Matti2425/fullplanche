function [T, P, rho, cp_mean] = compute_thermo(v_mag, NumericalParameters, Init, gas, S)
% COMPUTE_THERMO Computes thermodynamic properties from velocity field
%
% This function calculates the temperature, pressure, and density fields in a flow
% domain based on reference conditions and the velocity magnitude field.
% Uses isentropic flow relations and a polynomial model for specific heat capacity.
%
% INPUTS:
%   v_mag               - Velocity magnitude field at each grid point (N+2, M+2)
%   NumericalParameters - Structure containing N and M (number of nodes)
%   Init                - Structure containing reference T, P, rho, V
%   gas                 - Structure containing R (gas constant) and gamma (heat capacity ratio)
%   S                   - Logical array indicating solid regions
%
% OUTPUTS:
%   T       - Computed temperature field (N+2, M+2)
%   P       - Computed pressure field (N+2, M+2)
%   rho     - Computed density field (N+2, M+2)
%   cp_mean - Mean specific heat capacity field (N+2, M+2)
%
% NOTES:
%   - Uses polynomial model for cp: cp(T) = A + B*T + C*T^2
%   - Applies isentropic flow relations
%   - Fully vectorized implementation (no loops)
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    
    % Reference values
    T_ref = Init.T;
    P_ref = Init.P;
    v_ref = Init.V;
    R = gas.R;
    gamma = gas.gamma;
    
    % Polynomial coefficients for cp(T) = A + B*T + C*T^2
    A = 1005;      % J/(kg·K)
    B = 0.1;       % J/(kg·K^2)
    C = -0.0002;   % J/(kg·K^3)
    
    % Initialize arrays
    T = ones(N+2, M+2) * T_ref;
    P = ones(N+2, M+2) * P_ref;
    rho = ones(N+2, M+2) * Init.rho;
    
    % Iterative calculation for internal nodes (vectorized)
    % We need iteration because cp depends on T, and T depends on cp
    max_iter = 50;
    tol = 1e-6;
    
    for iter = 1:max_iter
        T_old = T;
        
        % Compute mean specific heat capacity (vectorized)
        % cp_mean = (1/(T-T_ref)) * integral(cp(T)dT from T_ref to T)
        dT = T - T_ref;
        
        % Avoid division by zero at reference temperature
        dT(abs(dT) < 1e-10) = 1e-10;
        
        cp_mean = (A .* dT + (B / 2) .* (T.^2 - T_ref^2) + (C / 3) .* (T.^3 - T_ref^3)) ./ dT;
        
        % Compute temperature from energy equation (vectorized)
        % T = T_ref + (v_ref^2 - v^2) / (2 * cp_mean)
        T = T_ref + (v_ref^2 - v_mag.^2) ./ (2 .* cp_mean);
        
        % Check convergence
        if max(abs(T(:) - T_old(:))) < tol
            break;
        end
    end
    
    % Compute pressure using isentropic relation (vectorized)
    % P = P_ref * (T / T_ref)^(gamma / (gamma - 1))
    P = P_ref .* (T ./ T_ref).^(gamma / (gamma - 1));
    
    % Compute density using ideal gas law (vectorized)
    % rho = P / (R * T)
    rho = P ./ (R .* T);
    
    % Apply boundary conditions for solid regions
    if nargin >= 5 && ~isempty(S)
        T(S) = T_ref;
        P(S) = P_ref;
        rho(S) = 0;  % Zero density in solid
        cp_mean(S) = A;  % Constant cp in solid
    end
    
end