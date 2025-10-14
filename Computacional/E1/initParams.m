function [psi, rho] = initParams(Init, NumericalParameters, solid)
% INITPARAMS Initializes the stream function and density fields
%
% This function initializes the computational fields with starting values
%
% Inputs:
%   Init                - Structure containing initial values (phi_start, rho)
%   NumericalParameters - Structure containing N and M
%   solid               - Matrix indicating solid (0) or fluid (1) regions
%
% Outputs:
%   psi - Initial stream function field (M+2, N+2)
%   rho - Initial density field (M+2, N+2)
%--------------------------------------------------------------------------

    N = NumericalParameters.N;
    M = NumericalParameters.M;
    
    % Initialize stream function
    psi = Init.phi_start * ones(M+2, N+2);
    
    % Initialize density (rho_F for fluid, 0 for solid)
    rho = Init.rho * solid;

end