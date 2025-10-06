function [v_x, v_y, v_mag] = compute_vel(phi, X, Y, NumericalParameters, Vin, Geometry, rho)
% COMPUTE_VEL Calculates the velocity field from the stream function
%
% This function calculates the velocity field in a computational domain based on the 
% stream function phi. It determines the velocity components in the x and y directions 
% and computes the velocity magnitude at each node.
%
% INPUTS:
%   phi                 - Stream function values at each computational node
%   X                   - X-coordinates of the grid points (N+2, M+2)
%   Y                   - Y-coordinates of the grid points (N+2, M+2)
%   NumericalParameters - Structure containing N and M (number of nodes)
%   Vin                 - Inlet velocity (free-stream velocity)
%   Geometry            - Structure containing rho_ref (reference density)
%   rho                 - Density field (N+2, M+2)
%
% OUTPUTS:
%   v_x   - Velocity component in the x-direction
%   v_y   - Velocity component in the y-direction
%   v_mag - Velocity magnitude at each node
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    rho_ref = Geometry.rho_ref;

    % Initialize velocity arrays
    v_x = zeros(N+2, M+2); 
    v_y = zeros(N+2, M+2);

    % Initialize density ratio arrays for blocking-off method
    rr_re = zeros(N+2, M+2);
    rr_rw = zeros(N+2, M+2);
    rr_rn = zeros(N+2, M+2);
    rr_rs = zeros(N+2, M+2);

    % Boundary conditions
    
    % Top and bottom walls → Tangential velocity = 0
    for i = 1:N+2
        v_y(i, 1) = 0;      % Bottom wall
        v_y(i, M+2) = 0;    % Top wall
    end

    % Developed flow at outlet
    for j = 1:M+2
        v_x(N+2, j) = v_x(N+1, j);
        v_y(N+2, j) = v_y(N+1, j);
    end
                
    % Channel inlet (x = 0) → Vin
    for j = 1:M+2
        v_x(1, j) = Vin;
        v_y(1, j) = 0;
    end

    % Compute velocities at internal nodes
    for i = 2:N+1
        for j = 2:M+1 
            
            % Blocking-off method: compute density ratios
            rr_re(i,j) = (X(i+1,j) - X(i,j)) / ...
                (((0.5 * (X(i+1,j) - X(i,j))) / (rho_ref / rho(i,j))) + ...
                 ((0.5 * (X(i+1,j) - X(i,j))) / (rho_ref / rho(i+1,j))));
            
            rr_rw(i,j) = (X(i,j) - X(i-1,j)) / ...
                (((0.5 * (X(i,j) - X(i-1,j))) / (rho_ref / rho(i,j))) + ...
                 ((0.5 * (X(i,j) - X(i-1,j))) / (rho_ref / rho(i-1,j))));
            
            rr_rn(i,j) = (Y(i,j+1) - Y(i,j)) / ...
                (((0.5 * (Y(i,j+1) - Y(i,j))) / (rho_ref / rho(i,j))) + ...
                 ((0.5 * (Y(i,j+1) - Y(i,j))) / (rho_ref / rho(i,j+1))));
            
            rr_rs(i,j) = (Y(i,j) - Y(i,j-1)) / ...
                (((0.5 * (Y(i,j) - Y(i,j-1))) / (rho_ref / rho(i,j))) + ...
                 ((0.5 * (Y(i,j) - Y(i,j-1))) / (rho_ref / rho(i,j-1))));
            
            % Velocity at faces (using stream function derivatives)
            % v_x = dphi/dy, v_y = -dphi/dx
            v_xn = +(1/rr_rn(i,j)) * (phi(i,j+1) - phi(i,j)) / (Y(i,j+1) - Y(i,j));
            v_xs = +(1/rr_rs(i,j)) * (phi(i,j) - phi(i,j-1)) / (Y(i,j) - Y(i,j-1));
            v_ye = -(1/rr_re(i,j)) * (phi(i+1,j) - phi(i,j)) / (X(i+1,j) - X(i,j));
            v_yw = -(1/rr_rw(i,j)) * (phi(i,j) - phi(i-1,j)) / (X(i,j) - X(i-1,j));
            
            % Velocity components at nodes (average of face velocities)
            v_x(i,j) = (v_xn + v_xs) / 2;
            v_y(i,j) = (v_ye + v_yw) / 2;

        end
    end
    
    % Velocity magnitude
    v_mag = sqrt(v_x.^2 + v_y.^2); 
    
end