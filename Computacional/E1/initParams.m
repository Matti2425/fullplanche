function [phi_start, T, P, rho] = initParams(Init, Numerical, S)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    phi_start = 0;
    M1 = ones(Numerical.N+2, Numerical.M+2);
    rho = ones(Numerical.N+2, Numerical.M+2)*Init.rho;
    T = M1*Init.T;
    P = M1*Init.P;

    rho(S) = 0;

end