function [Geometry,Vin,Pin,Tin,phi_body,gas,Numerical,Init] = readInputData()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    Geometry.L = 20;
    Geometry.H = 5;
    Geometry.D = 1;
    Geometry.Center = [Geometry.L/4,Geometry.H/2];
    Vin = 1;
    Pin = 101325;
    Tin = 1;
    phi_body = Vin/2;
    gas.cp = 1;
    gas.gamma = 1.4;
    gas.R = 281;
    Numerical.N = 10;
    Numerical.M = 5;
    Geometry.rho_ref = 1.225;

    Init.T = Tin;
    Init.P = Pin;
    Init.rho = Geometry.rho_ref;
    Init.V = Vin;
end