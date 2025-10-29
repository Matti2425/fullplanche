function [a_n, a_s, a_e, a_w, a_p, a_p_base, b_p] = computeCoefficients(v_x, v_y, Numerical, Fluid, ...
                                                                      sx_w, sx_e, sx_p, sy_s, sy_n, sy_p, ...
                                                                      gamma_val, case_study, A_x, A_y)
% COMPUTECOEFFICIENTS Evaluates discretization coefficients for convection-diffusion equation
%
% This function computes the coefficients for the finite volume discretization
% of the generic convection-diffusion transport equation using various schemes.
%
% INPUTS:
%   v_x, v_y   - Velocity field components
%   Numerical  - Structure with N, M, time_step
%   Fluid      - Structure with rho
%   sx_w, sx_e, sx_p - Geometric parameters in x-direction
%   sy_s, sy_n, sy_p - Geometric parameters in y-direction
%   gamma_val  - Current diffusion coefficient value
%   case_study - Discretization scheme (1=Upwind, 2=Exp, 3=CDS, 4=HDS, 5=PLDS)
%   A_x, A_y   - Cell face areas
%
% OUTPUTS:
%   a_n, a_s, a_e, a_w - Neighbor coefficients (north, south, east, west)
%   a_p                - Central coefficient
%   a_p_base           - Transient term coefficient
%   b_p                - Source term
%--------------------------------------------------------------------------

    N = Numerical.N;
    M = Numerical.M;
    rho = Fluid.rho;
    time_step = Numerical.time_step;
    
    % Initialize coefficient matrices (N+2, M+2)
    a_n = zeros(N+2, M+2);
    a_s = zeros(N+2, M+2);
    a_e = zeros(N+2, M+2);
    a_w = zeros(N+2, M+2);
    a_p = ones(N+2, M+2);
    a_p_base = zeros(N+2, M+2);
    b_p = zeros(N+2, M+2);
    
    % Compute diffusion coefficients (vectorized where possible)
    D_n = (gamma_val * A_x) ./ sy_n(2:M+1);
    D_s = (gamma_val * A_x) ./ sy_s(2:M+1);
    D_e = (gamma_val * A_y) ./ sx_e(2:N+1);
    D_w = (gamma_val * A_y) ./ sx_w(2:N+1);
    
    % Loop over interior cells
    for j = 2:M+1
        for i = 2:N+1
            % Compute convective fluxes at cell faces
            F_n = rho * (v_y(i, j+1) + v_y(i, j)) * A_x / 2;
            F_s = rho * (v_y(i, j-1) + v_y(i, j)) * A_x / 2;
            F_e = rho * (v_x(i+1, j) + v_x(i, j)) * A_y / 2;
            F_w = rho * (v_x(i-1, j) + v_x(i, j)) * A_y / 2;
            
            % Compute Peclet numbers
            p_n = F_n / D_n(j-1);
            p_s = F_s / D_s(j-1);
            p_e = F_e / D_e(i-1);
            p_w = F_w / D_w(i-1);
            
            % Select discretization scheme
            switch case_study
                case "UDS" 
                    A_n = 1;
                    A_s = 1;
                    A_e = 1;
                    A_w = 1;

                case "CDS"  
                    A_n = abs((abs(p_n)/2) - 1);
                    A_s = abs((abs(p_s)/2) - 1);
                    A_e = abs((abs(p_e)/2) - 1);
                    A_w = abs((abs(p_w)/2) - 1);

                case "HDS" 
                    A_n = max(0, 1 - (abs(p_n)/2));
                    A_s = max(0, 1 - (abs(p_s)/2));
                    A_e = max(0, 1 - (abs(p_e)/2));
                    A_w = max(0, 1 - (abs(p_w)/2));
                    
                case "EDS" 
                    A_n = abs(p_n) / (exp(abs(p_n)) - 1);
                    A_s = abs(p_s) / (exp(abs(p_s)) - 1);
                    A_e = abs(p_e) / (exp(abs(p_e)) - 1);
                    A_w = abs(p_w) / (exp(abs(p_w)) - 1);
                    
                case "PLDS" 
                    A_n = max(0, 1 - (abs(p_n)/2)^5);
                    A_s = max(0, 1 - (abs(p_s)/2)^5);
                    A_e = max(0, 1 - (abs(p_e)/2)^5);
                    A_w = max(0, 1 - (abs(p_w)/2)^5);
                    
                otherwise
                    error('The case_study selected is not available');
            end
            
            % Compute coefficients using the selected scheme
            a_n(i, j) = D_n(j-1) * A_n + max(-F_n, 0);
            a_s(i, j) = D_s(j-1) * A_s + max(F_s, 0);
            a_e(i, j) = D_e(i-1) * A_e + max(-F_e, 0);
            a_w(i, j) = D_w(i-1) * A_w + max(F_w, 0);
            
            % Transient term
            a_p_base(i, j) = rho * sx_p(i) * sy_p(j) / time_step;
            
            % Central coefficient
            a_p(i, j) = a_e(i, j) + a_w(i, j) + a_n(i, j) + a_s(i, j) + a_p_base(i, j);
        end
    end

end