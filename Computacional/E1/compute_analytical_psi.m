function psi_analytical = compute_analytical_psi(x_P, y_P, Geometry, v_in, CASE, w, R_cyl)
% COMPUTE_ANALYTICAL_PSI Computes the analytical stream function based on the case
%
% This function computes the analytical stream function for three cases:
% - CASE 1: Uniform channel flow (Eq 5.1.1)
% - CASE 2: Flow around static cylinder with source/sink (Eq 5.2.1)
% - CASE 3: Flow around rotating cylinder with circulation (Eq 5.3.1)
%
% INPUTS:
%   x_P, y_P        - 1D arrays of mesh coordinates
%   Geometry        - Structure with ObjectCenter, H, L, D
%   v_in            - Inlet velocity
%   CASE            - Case selector (1, 2, or 3)
%   w               - Angular velocity of cylinder (for CASE 3)
%   R_cyl           - Cylinder radius
%
% OUTPUTS:
%   psi_analytical  - Analytical stream function field (M+2, N+2)
%--------------------------------------------------------------------------

    M = length(y_P) - 2;
    N = length(x_P) - 2;
    psi_analytical = zeros(M+2, N+2);
    
    x0 = Geometry.ObjectCenter(1);
    y0 = Geometry.ObjectCenter(2);
    H = Geometry.H;
    L = Geometry.L;
    
    % Create meshgrid for computation
    [X, Y] = meshgrid(x_P, y_P);
    
    switch CASE
        case 1
            % Equation 5.1.1: Uniform channel flow
            % psi(x, y) = -V*x + U*y
            % For uniform flow in x-direction: set U = v_in, V = 0
            % This gives psi = v_in*y which matches inlet BC
            psi_analytical = v_in * Y;
            
        case 2
            % Equation 5.2.1: Flow around static cylinder (source/doublet)
            % psi(x, y) = -V*x + U*y + (mu/2*pi) * ((y-y0)*cos(alpha) - (x-x0)*sin(alpha)) / ((x-x0)^2 + (y-y0)^2)
            % 
            % Physical interpretation:
            % - Uniform flow: psi = v_in * y
            % - Perturbation from cylinder: doublet term
            % - Alpha = 0: doublet aligned with horizontal axis
            
            % Compute distance from cylinder center
            r_sq = (X - x0).^2 + (Y - y0).^2;
            r = sqrt(r_sq);
            
            % Uniform flow part
            psi_uniform = v_in * Y;
            
            % For potential flow around a cylinder: doublet strength
            % mu = 2*pi*R_cyl^2*v_in represents the strength of the doublet
            mu = 2 * pi * R_cyl^2 * v_in;
            
            % Doublet part (with alpha = 0: cos(0)=1, sin(0)=0)
            % psi_doublet = (mu/2*pi) * (Y - y0) / r_sq
            psi_doublet = zeros(size(X));
            psi_doublet(r > R_cyl) = (mu / (2*pi)) * (Y(r > R_cyl) - y0) ./ r_sq(r > R_cyl);
            
            % Total analytical solution
            psi_analytical = psi_uniform + psi_doublet;
            
            % Inside cylinder: psi = constant (psi_body)
            psi_analytical(r <= R_cyl) = v_in * y0;
            
        case 3
            % Equation 5.3.1: Flow around rotating cylinder (with circulation)
            % psi(x, y) = psi_static - (Gamma/2*pi) * ln(sqrt((x-x0)^2 + (y-y0)^2))
            % 
            % psi_static includes: uniform flow + doublet from cylinder
            % plus circulation term for rotating cylinder
            
            % Compute distance from cylinder center
            r = sqrt((X - x0).^2 + (Y - y0).^2);
            r_sq = r.^2;
            
            % Static component (uniform + doublet)
            mu = 2 * pi * R_cyl^2 * v_in;
            psi_uniform = v_in * Y;
            
            psi_static = psi_uniform;
            psi_static(r > R_cyl) = psi_static(r > R_cyl) + (mu / (2*pi)) * (Y(r > R_cyl) - y0) ./ r_sq(r > R_cyl);
            psi_static(r <= R_cyl) = v_in * y0;
            
            % Circulation component: Gamma related to angular velocity
            % Gamma = 2*pi*R_cyl^2*w represents vortex strength
            Gamma = 2 * pi * R_cyl^2 * w;
            psi_circulation = zeros(size(X));
            psi_circulation(r > R_cyl) = -(Gamma / (2*pi)) * log(r(r > R_cyl) / R_cyl);
            psi_circulation(r <= R_cyl) = 0;  % Inside cylinder, circulation is not added
            
            % Total analytical solution
            psi_analytical = psi_static + psi_circulation;
    end
    
end