function [phi, T, P, rho] = initParams(Init, Numerical, S, phi_body)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    phi = ones(Numerical.N+2, Numerical.M+2);
    M1 = ones(Numerical.N+2, Numerical.M+2);
    T = M1*Init.T;
    P = M1*Init.P;
    rho = M1*Init.rho;

    phi(logical(S)) = phi_body;


end