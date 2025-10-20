function [domain, inlet, boundary, material, solver, state, rotSpeed] = setupSimulation(simCase)
% SETUPSIMULATION Initializes all simulation parameters for potential flow analysis
%
% This function configures all physical, geometrical, and numerical parameters
% necessary for solving potential flow around a cylindrical obstacle.
%
% INPUTS:
%   simCase   - Selected Case
%
% OUTPUTS:
%   domain    - Geometric domain structure (length, height, diameter, center, ref_density)
%   inlet     - Inlet conditions structure (velocity)
%   boundary  - Boundary conditions structure (psi_surface)
%   material  - Fluid material properties (gamma, specific_gas_constant)
%   solver    - Numerical solver parameters (nx, ny, tolerance, max_iterations)
%   state     - Initial state structure (temperature, pressure, density, velocity, psi_ref)
%   rotSpeed  - Angular velocity of rotating cylinder (rad/s)
%   simCase   - Case selector (1=channel, 2=static, 3=rotating)
%--------------------------------------------------------------------------

    % Domain Geometry Configuration
    domain.length = 2;              % Channel length [m]
    domain.height = 1;              % Channel height [m]
    domain.diameter = 0.2;          % Obstacle diameter [m]
    domain.center_x = domain.length / 2;  % Obstacle x-position [m]
    domain.center_y = domain.height / 2;  % Obstacle y-position [m]

    % Fluid Material Properties
    material.gamma = 1.4;           % Heat capacity ratio (air)
    material.R_gas = 287;  % Gas constant [J/(kg·K)]

    % Flow Inlet Conditions
    inlet.velocity = 1;             % Inlet velocity [m/s]
    inlet.pressure = 101325;        % Reference pressure [Pa]
    inlet.temperature = 300;        % Reference temperature [K]
    
    % Compute reference density
    rho_reference = inlet.pressure / (material.R_gas * inlet.temperature);
    domain.ref_density = rho_reference;

    % Cylinder Rotation Parameter
    if simCase == 3
        rotSpeed = 25;  % Angular velocity for rotating cylinder [rad/s]
    else
        rotSpeed = 0;   % Stationary cylinder
    end

    % Boundary Condition on Obstacle Surface
    if simCase ~= 1
        R_obstacle = domain.diameter / 2;
        boundary.psi_surface = inlet.velocity * (domain.height/2) + ...
                              R_obstacle^2 * rotSpeed;
    else
        boundary.psi_surface = 0;
    end

    % Numerical Solver Configuration
    solver.nx = 100;                % Grid points x-direction
    solver.ny = 100;                % Grid points y-direction
    solver.tolerance = 1e-10;       % Convergence tolerance
    solver.max_iterations = 1000000;  % Maximum iteration count

    % Initial Conditions
    state.temperature = inlet.temperature;
    state.pressure = inlet.pressure;
    state.density = rho_reference;
    state.velocity = inlet.velocity;
    state.psi_ref = 0;              % Reference stream function value

end