function [psi, density] = initializeFields(state, solver, fluidMask, varargin)
% INITIALIZEFIELDS Sets initial conditions for stream function and density
%
% Initializes computational fields with reference values and applies
% boundary conditions for solid/fluid regions.
%
% Inputs:
%   state     - Structure with initial values (psi_ref, density)
%   solver    - Structure containing nx and ny
%   fluidMask - Binary matrix indicating fluid (1) or solid (0) regions
%
% Outputs:
%   streamFunc - Initial stream function field (ny+2, nx+2)
%   density    - Initial density field (ny+2, nx+2)
%--------------------------------------------------------------------------

    nx = solver.nx;
    ny = solver.ny;

    % Initialize stream function with reference value everywhere
    psi = state.psi_ref * ones(ny+2, nx+2);

    % Initialize density: reference density in fluid, zero in solid
    density = state.density * fluidMask;

end