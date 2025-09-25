function [phi, T, P, rho] = initParams(Init, Numerical, S, phi_body)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    phi = ones(Numerical.N, Numerical.M);
    M1 = ones(Numerical.N, Numerical.M);
    T = M1*Init.T;
    P = M1*Init.P;
    rho = M1*Init.rho;

    phi(S) = phi_body;


end