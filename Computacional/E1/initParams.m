function [phi, T, P, rho] = initParams(Init, Numerical, S)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    phi = ones(Numerical.N+2,Numerical.M+2)*Init.phi_start;
    M1 = ones(Numerical.N+2, Numerical.M+2);
    rho = ones(Numerical.N+2, Numerical.M+2)*Init.rho;
    T = M1*Init.T;
    P = M1*Init.P;

    rho(S) = 0;

end