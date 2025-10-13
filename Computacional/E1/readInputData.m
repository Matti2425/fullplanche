function [Geometry,Vin,Pin,Tin,phi_body,gas,Numerical,Init] = readInputData()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Geometry properties
    Geometry.L = 2;
    Geometry.H = 1;
    Geometry.D = 0.2;
    Geometry.ObjectCenter = [Geometry.L/2,Geometry.H/2];

    % Fluid properties
    gas.cp = 1;
    gas.gamma = 1.4;
    gas.R = 287;
    Geometry.rho_ref = 1.225;
        
    % Input values
    Vin = 1;
    Pin = 101325;
    Tin = 300;
    Rho_in = Pin/(gas.R*Tin);

    % Body properties
    phi_body = Vin*Geometry.H/2;


    % Numerical system properties
    Numerical.N = 5;
    Numerical.M = 7;
    Numerical.d = 1e-7;
    Numerical.i_max = 1000;

    % Initial properties
    Init.T = Tin;
    Init.P = Pin;
    Init.rho = Rho_in;
    Init.V = Vin;
    Init.phi_start = 0;

end