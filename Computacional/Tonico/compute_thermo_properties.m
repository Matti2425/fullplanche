function [T_P, p_P, rho_P, cp_mean] = compute_thermo_properties(T_ref, p_ref, v_ref, v_P, R, gamma, M, N)

% COMPUTE_THERMO_PROPERTIES Computes thermodynamic properties from velocity and reference values.
%
% This function calculates the temperature, pressure, and density fields in a flow
% domain based on a reference temperature, pressure, and velocity.
%
% Inputs:
%   T_ref   - Reference temperature.
%   p_ref   - Reference pressure.
%   v_ref   - Reference velocity.
%   v_P     - Velocity magnitude field at each grid point.
%   R       - Specific gas constant.
%   gamma   - Heat capacity ratio (Cp/Cv).
%   M, N    - Number of grid points in the y and x directions.
%
% Outputs:
%   T_P     - Computed temperature field.
%   p_P     - Computed pressure field.
%   rho_P   - Computed density field.
%   cp_mean - Computed specific heat capacity field.
%
% A polynomial model is used to estimate the specific heat capacity as a function of temperature.
%--------------------------------------------------------------------------

A = 1005;     
B = 0.1;       
C = -0.0002;  

T_P = zeros(M+2, N+2);
p_P = zeros(M+2, N+2);
rho_P = zeros(M+2, N+2);
 
cp_mean = (1 ./ (T_P - T_ref)) .*(A .* (T_P - T_ref) +(B / 2) .* (T_P.^2 - T_ref^2) +(C / 3) .* (T_P.^3 - T_ref^3));
    
for i = 2:N+1
    for j = 2:M+1
        T_P(j,i) = T_ref + ((v_ref^2 - v_P(j,i)^2) / (2 * cp_mean(j,i)));
        p_P(j,i) = p_ref * (T_P(j,i) / T_ref)^(gamma / (gamma - 1));
        rho_P(j,i) = p_P(j,i) / (R * T_P(j,i));
        end
    end
end
