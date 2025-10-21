function [T_P, p_P, rho_P, cp_mean] = compute_thermo(v_P, NumericalParameters, Init, gas)
% COMPUTE_THERMO Computes thermodynamic properties from velocity field
%
% This function calculates the temperature, pressure, and density fields in a flow
% domain based on reference conditions and the velocity magnitude field.
% Uses isentropic flow relations and a polynomial model for specific heat capacity.
%
% INPUTS:
%   v_P                 - Velocity magnitude field at each grid point (M+2, N+2)
%   NumericalParameters - Structure containing N and M (number of nodes)
%   Init                - Structure containing reference T, P, rho, V
%   gas                 - Structure containing R (gas constant) and gamma (heat capacity ratio)
%
% OUTPUTS:
%   T_P     - Computed temperature field (M+2, N+2)
%   p_P     - Computed pressure field (M+2, N+2)
%   rho_P   - Computed density field (M+2, N+2)
%   cp_mean - Mean specific heat capacity field (M+2, N+2)
%
% NOTES:
%   - Uses polynomial model for cp: cp(T) = A + B*T + C*T^2
%   - Applies isentropic flow relations
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    
    % Reference values
    T_ref = Init.T;
    p_ref = Init.P;
    v_ref = Init.V;
    R = gas.R;
    gamma = gas.gamma;
    
    % Polynomial coefficients for cp(T) = A + B*T + C*T^2
    A = 1005;      % J/(kg·K)
    B = 0.1;       % J/(kg·K^2)
    C = -0.0002;   % J/(kg·K^3)
    
    % Initialize arrays (M+2, N+2)
    T_P = zeros(M+2, N+2);
    p_P = zeros(M+2, N+2);
    rho_P = zeros(M+2, N+2);
    
    % Compute mean specific heat capacity (vectorized)
    % cp_mean = (1/(T-T_ref)) * integral(cp(T)dT from T_ref to T)
    cp_mean = (1 ./ (T_P - T_ref)) .* (A .* (T_P - T_ref) + ...
              (B / 2) .* (T_P.^2 - T_ref^2) + ...
              (C / 3) .* (T_P.^3 - T_ref^3));
    
    % Compute thermodynamic properties for internal nodes
    for i = 2:N+1
        for j = 2:M+1
            T_P(j,i) = T_ref + ((v_ref^2 - v_P(j,i)^2) / (2 * cp_mean(j,i)));
            p_P(j,i) = p_ref * (T_P(j,i) / T_ref)^(gamma / (gamma - 1));
            rho_P(j,i) = p_P(j,i) / (R * T_P(j,i));
        end
    end
    
end