function [tempField, pressField, rhoField, cpAvg] = thermoProperties(velMag, solver, state, material)
% THERMOPROPERTIES Calculates thermodynamic properties from flow field
%
% Computes temperature, pressure, and density using isentropic flow relations
% and polynomial model for specific heat capacity.
%
% INPUTS:
%   velMag  - Velocity magnitude field at grid points (ny+2, nx+2)
%   solver  - Structure containing nx and ny
%   state   - Structure with reference T, P, rho, V
%   material - Structure with R (gas constant) and gamma (heat ratio)
%
% OUTPUTS:
%   tempField - Temperature field (ny+2, nx+2)
%   pressField - Pressure field (ny+2, nx+2)
%   rhoField   - Density field (ny+2, nx+2)
%   cpAvg      - Average specific heat capacity field (ny+2, nx+2)
%--------------------------------------------------------------------------

    nx = solver.nx;
    ny = solver.ny;

    % Reference state values
    tempRef = state.temperature;
    pressRef = state.pressure;
    velRef = state.velocity;
    R_gas = material.specific_gas_constant;
    gammaRatio = material.gamma;

    % Polynomial coefficients for cp(T) = A + B*T + C*T^2
    cp_A = 1005;        % J/(kg·K) - intercept
    cp_B = 0.1;         % J/(kg·K^2) - linear coefficient
    cp_C = -0.0002;     % J/(kg·K^3) - quadratic coefficient

    % Initialize arrays
    tempField = zeros(ny+2, nx+2);
    pressField = zeros(ny+2, nx+2);
    rhoField = zeros(ny+2, nx+2);
    cpAvg = zeros(ny+2, nx+2);

    % Compute thermodynamic properties at interior nodes
    for i = 2:nx+1
        for j = 2:ny+1
            % Temperature using energy equation
            % T = T_ref + (V_ref^2 - V^2) / (2 * cp_mean)
            % First estimate cp_mean assuming small temperature change
            cp_local = cp_A + cp_B * tempRef + cp_C * tempRef^2;

            % Temperature from isentropic relation
            tempField(j,i) = tempRef + ((velRef^2 - velMag(j,i)^2) / (2 * cp_local));

            % Refine cp_mean using computed temperature
            tempMid = (tempField(j,i) + tempRef) / 2;
            cpAvg(j,i) = (1 / (tempField(j,i) - tempRef)) * ...
                        (cp_A * (tempField(j,i) - tempRef) + ...
                         (cp_B / 2) * (tempField(j,i)^2 - tempRef^2) + ...
                         (cp_C / 3) * (tempField(j,i)^3 - tempRef^3));

            % Recompute temperature with refined cp (iterative refinement)
            tempField(j,i) = tempRef + ((velRef^2 - velMag(j,i)^2) / (2 * cpAvg(j,i)));

            % Pressure from isentropic relation
            pressField(j,i) = pressRef * (tempField(j,i) / tempRef)^(gammaRatio / (gammaRatio - 1));

            % Density from ideal gas law
            rhoField(j,i) = pressField(j,i) / (R_gas * tempField(j,i));
        end
    end

    % Apply boundary conditions to scalar fields (copy interior to boundaries)
    % This ensures physical consistency at domain edges
    tempField(1,:) = tempRef;
    tempField(ny+2,:) = tempRef;
    tempField(:,1) = tempRef;
    tempField(:,nx+2) = tempRef;

    pressField(1,:) = pressRef;
    pressField(ny+2,:) = pressRef;
    pressField(:,1) = pressRef;
    pressField(:,nx+2) = pressRef;

    rhoField(1,:) = state.density;
    rhoField(ny+2,:) = state.density;
    rhoField(:,1) = state.density;
    rhoField(:,nx+2) = state.density;

end