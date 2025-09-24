function [Geometry,vin,pin,Tin,phi_body,gas,Numerical] = readInputData()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    Geometry.L = 20;
    Geometry.H = 5;
    Geometry.D = 1;
    Geometry.Center = [Geometry.L/4,Geometry.H/2];
    vin = 1;
    pin = 101325;
    Tin = 1;
    phi_body = vin/2;
    gas.cp = 1;
    gas.gamma = 1.4;
    gas.R = 281;
    Numerical.N = 600;
    Numerical.M = 200;
end