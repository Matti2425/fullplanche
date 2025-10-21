function [velX, velY, velMag] = velocityField(streamFunc, gridX, gridY, solver, velocity, rho_ref, rho)
% VELOCITYFIELD Reconstructs velocity components from stream function
%
% Computes velocity field using the stream function through finite
% difference derivatives, applying the blocking-off method for density variation.
%
% INPUTS:
%   streamFunc  - Stream function values at grid points (ny+2, nx+2)
%   gridX, gridY - Coordinate vectors
%   solver      - Structure containing nx and ny
%   velocity    - Inlet velocity
%   rho_ref     - Reference density
%   rho         - Density field (ny+2, nx+2)
%
% OUTPUTS:
%   velX   - X-component of velocity (ny+2, nx+2)
%   velY   - Y-component of velocity (ny+2, nx+2)
%   velMag - Magnitude of velocity (ny+2, nx+2)
%--------------------------------------------------------------------------

    nx = solver.nx;
    ny = solver.ny;

    % Initialize velocity arrays
    velX = zeros(ny+2, nx+2);
    velY = zeros(ny+2, nx+2);

    % Precompute density ratios for blocking-off method
    rho_ratio_east = zeros(ny+2, nx+2);
    rho_ratio_west = zeros(ny+2, nx+2);
    rho_ratio_north = zeros(ny+2, nx+2);
    rho_ratio_south = zeros(ny+2, nx+2);

    % Boundary conditions: No normal velocity at walls
    for i = 1:nx+2
        velY(1, i) = 0;      % Bottom wall
        velY(ny+2, i) = 0;   % Top wall
    end

    % Outlet condition: Developed flow (zero gradient)
    for j = 1:ny+2
        velX(j, nx+2) = velX(j, nx+1);
        velY(j, nx+2) = velY(j, nx+1);
    end

    % Inlet condition: Uniform velocity
    for j = 1:ny+2
        velX(j, 1) = velocity;
        velY(j, 1) = 0;
    end

    % Compute velocity at interior nodes
    for i = 2:nx+1
        for j = 2:ny+1

            % Compute density ratios (blocking-off method)
            dx_east = gridX(i+1) - gridX(i);
            dx_west = gridX(i) - gridX(i-1);
            dy_north = gridY(j+1) - gridY(j);
            dy_south = gridY(j) - gridY(j-1);

            rho_ratio_east(j,i) = dx_east / ...
                (((0.5 * dx_east) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * dx_east) / (rho_ref / rho(j,i+1))));

            rho_ratio_west(j,i) = dx_west / ...
                (((0.5 * dx_west) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * dx_west) / (rho_ref / rho(j,i-1))));

            rho_ratio_north(j,i) = dy_north / ...
                (((0.5 * dy_north) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * dy_north) / (rho_ref / rho(j+1,i))));

            rho_ratio_south(j,i) = dy_south / ...
                (((0.5 * dy_south) / (rho_ref / rho(j,i))) + ...
                 ((0.5 * dy_south) / (rho_ref / rho(j-1,i))));

            % Compute velocity at cell faces
            % v_x = dpsi/dy, v_y = -dpsi/dx
            vel_north = +(1/rho_ratio_north(j,i)) * (streamFunc(j+1, i) - streamFunc(j, i)) / dy_north;
            vel_south = +(1/rho_ratio_south(j,i)) * (streamFunc(j, i) - streamFunc(j-1, i)) / dy_south;
            vel_east = -(1/rho_ratio_east(j,i)) * (streamFunc(j, i+1) - streamFunc(j, i)) / dx_east;
            vel_west = -(1/rho_ratio_west(j,i)) * (streamFunc(j, i) - streamFunc(j, i-1)) / dx_west;

            % Average face velocities to get node-centered values
            velX(j, i) =  (vel_north + vel_south) / 2;  
            velY(j, i) = (vel_east + vel_west) / 2;     

        end
    end

    % Compute velocity magnitude
    velMag = sqrt(velX.^2 + velY.^2);

end