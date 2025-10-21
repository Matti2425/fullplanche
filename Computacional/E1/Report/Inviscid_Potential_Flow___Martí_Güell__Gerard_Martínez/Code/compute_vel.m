function [v_xP, v_yP, v_P] = compute_vel(psi, x_P, y_P, NumericalParameters, v_in, rho_ref, rho)
% COMPUTE_VEL Calculates the velocity field from the stream function
%
% This function calculates the velocity field in a computational domain based on the 
% stream function psi. It determines the velocity components in the x and y directions 
% and computes the velocity magnitude at each node.
%
% INPUTS:
%   psi                 - Stream function values at each computational node (M+2, N+2)
%   x_P, y_P            - 1D arrays of grid coordinates
%   NumericalParameters - Structure containing N and M (number of nodes)
%   v_in                - Inlet velocity (free-stream velocity)
%   rho_ref             - Reference density
%   rho                 - Density field (M+2, N+2)
%
% OUTPUTS:
%   v_xP  - Velocity component in the x-direction (M+2, N+2)
%   v_yP  - Velocity component in the y-direction (M+2, N+2)
%   v_P   - Velocity magnitude at each node (M+2, N+2)
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;

    % Initialize velocity arrays (M+2, N+2)
    v_xP = zeros(M+2, N+2); 
    v_yP = zeros(M+2, N+2);

    % Initialize density ratio arrays for blocking-off method
    rr_re = zeros(M+2, N+2);
    rr_rw = zeros(M+2, N+2);
    rr_rn = zeros(M+2, N+2);
    rr_rs = zeros(M+2, N+2);

    % Boundary conditions
    
    % Top and bottom walls → Tangential velocity = 0
    for i = 1:N+2
        v_yP(1, i) = 0;      % Bottom wall
        v_yP(M+2, i) = 0;    % Top wall
    end

    % Developed flow at outlet
    for j = 1:M+2
        v_xP(j, N+2) = v_xP(j, N+1);
        v_yP(j, N+2) = v_yP(j, N+1);
    end
                
    % Channel inlet (x = 0) → v_in
    for j = 1:M+2
        v_xP(j, 1) = v_in;
        v_yP(j, 1) = 0;
    end

    % Compute velocities at internal nodes
    for i = 2:N+1
        for j = 2:M+1 
            
            % Blocking-off method: compute density ratios
            rr_re(j,i) = (x_P(i+1) - x_P(i)) / ...
                (((0.5 * (x_P(i+1) - x_P(i))) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * (x_P(i+1) - x_P(i))) / (rho_ref / rho(j,i+1))));
            
            rr_rw(j,i) = (x_P(i) - x_P(i-1)) / ...
                (((0.5 * (x_P(i) - x_P(i-1))) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * (x_P(i) - x_P(i-1))) / (rho_ref / rho(j,i-1))));
            
            rr_rn(j,i) = (y_P(j+1) - y_P(j)) / ...
                (((0.5 * (y_P(j+1) - y_P(j))) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * (y_P(j+1) - y_P(j))) / (rho_ref / rho(j+1,i))));
            
            rr_rs(j,i) = (y_P(j) - y_P(j-1)) / ...
                (((0.5 * (y_P(j) - y_P(j-1))) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * (y_P(j) - y_P(j-1))) / (rho_ref / rho(j-1,i))));
            
            % Velocity at faces (using stream function derivatives)
            % v_x = dpsi/dy, v_y = -dpsi/dx
            v_xn = +(1/rr_rn(j,i)) * (psi(j+1, i) - psi(j, i)) / (y_P(j+1) - y_P(j));
            v_xs = +(1/rr_rs(j,i)) * (psi(j, i) - psi(j-1, i)) / (y_P(j) - y_P(j-1));
            v_ye = -(1/rr_re(j,i)) * (psi(j, i+1) - psi(j, i)) / (x_P(i+1) - x_P(i));
            v_yw = -(1/rr_rw(j,i)) * (psi(j, i) - psi(j, i-1)) / (x_P(i) - x_P(i-1));
            
            % Velocity components at nodes (average of face velocities)
            v_xP(j, i) = (v_xn + v_xs) / 2;
            v_yP(j, i) = (v_ye + v_yw) / 2;

        end
    end
    
    % Velocity magnitude (vectorized)
    v_P = sqrt(v_xP.^2 + v_yP.^2); 
    
end