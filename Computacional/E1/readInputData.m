function [Geometry,Vin,Pin,Tin,phi_body,gas,Numerical,Init] = readInputData()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    Geometry.L = 20;
    Geometry.H = 5;
    Geometry.D = 1;
    Geometry.ObjectCenter = [Geometry.L/4,Geometry.H/2];
    Vin = 1;
    Pin = 101325;
    Tin = 1;
    phi_body = Vin*Geometry.H/2;
    gas.cp = 1;
    gas.gamma = 1.4;
    gas.R = 281;
    Rho_in = Pin/(gas.R*Tin);
    Numerical.N = 10;
    Numerical.M = 5;
    Numerical.d = 1e-4;
    Geometry.rho_ref = 1.225;

    Init.T = Tin;
    Init.P = Pin;
    Init.rho = Rho_in;
    Init.V = Vin;
end