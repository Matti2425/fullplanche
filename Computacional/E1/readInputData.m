function [Geometry, v_in, p_in, T_in, psi_body, gas, Numerical, Init, w, CASE] = readInputData()
% READINPUTDATA Reads and sets up all input parameters for the simulation
%
% This function defines all the physical, geometrical, and numerical parameters
% required for the potential flow simulation around a cylinder.
%
% OUTPUTS:
%   Geometry - Structure with L, H, D, ObjectCenter, rho_ref
%   v_in     - Inlet velocity
%   p_in     - Inlet pressure
%   T_in     - Inlet temperature
%   psi_body - Stream function value on cylinder surface
%   gas      - Structure with R, gamma
%   Numerical - Structure with N, M, d (tolerance), i_max
%   Init     - Structure with initial values
%   w        - Angular velocity of cylinder (rad/s)
%   CASE     - Case selector (1=channel, 2=static cylinder, 3=rotating cylinder)
%--------------------------------------------------------------------------

    % CASE SELECTION
    % CASE=1 for channel only
    % CASE=2 for channel + static cylinder
    % CASE=3 for channel + rotating cylinder (w>0 --> Clockwise)
    CASE = 2;
    
    % Geometry properties
    Geometry.L = 2;
    Geometry.H = 1;
    Geometry.D = 0.2;
    Geometry.ObjectCenter = [Geometry.L/2, Geometry.H/2];

    % Fluid properties
    gas.gamma = 1.4;
    gas.R = 287;
        
    % Input values
    v_in = 1;
    p_in = 101325;
    T_in = 300;
    rho_in = p_in / (gas.R * T_in);
    Geometry.rho_ref = rho_in;

    % Cylinder angular velocity
    if CASE == 3
        w = 25;  % Angular velocity for rotating cylinder
    else
        w = 0;   % Static cylinder or no cylinder
    end

    % Body properties (stream function on cylinder surface)
    if CASE ~= 1
        R_cyl = Geometry.D / 2;
        psi_body = v_in * (Geometry.H/2) + R_cyl^2 * w * (1/(2*pi));
    else
        psi_body = 0;
    end

    % Numerical system properties
    Numerical.N = 50;
    Numerical.M = 50;
    Numerical.d = 1e-10;
    Numerical.i_max = 1000000;

    % Initial properties
    Init.T = T_in;
    Init.P = p_in;
    Init.rho = rho_in;
    Init.V = v_in;
    Init.phi_start = 0;

end