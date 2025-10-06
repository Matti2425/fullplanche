function [v_xP, v_yP, v_P] = compute_velocity(psi, y_P, x_P, M, N, v_in, rho_ref, rho)

% This function calculates the velocity field in a computational domain based on the 
% stream function ψ. It determines the velocity components in the x and y directions 
% and computes the velocity magnitude at each node.
%
% INPUTS:
%   psi     - Stream function values at each computational node.
%   y_P     - Y-coordinates of the grid points.
%   x_P     - X-coordinates of the grid points.
%   M, N    - Number of nodes in the Y and X directions, respectively.
%   v_in    - Inlet velocity (free-stream velocity).
%   rho_ref - Reference density.
%   rho     - Density field
%
% OUTPUTS:
%   v_xP - Velocity component in the x-direction.
%   v_yP - Velocity component in the y-direction.
%   v_P  - Velocity magnitude at each node.
%--------------------------------------------------------------------------

    v_xP = zeros(M+2, N+2); 
    v_yP = zeros(M+2, N+2);

    rr_re = zeros(M+2, N+2);
    rr_rw = zeros(M+2, N+2);
    rr_rn = zeros(M+2, N+2);
    rr_rs = zeros(M+2, N+2);


% top and bottom walls → Tangencial velocity=0
for i = 1:N+2
    v_yP(1, i) = 0;
    v_yP(M+2, i) = 0;

end

% Developed flow at output
for j = 1:M+2
    v_xP(j, N+2) = v_xP(j, N+1);
    v_yP(j, N+2) = v_yP(j, N+1);
end
            
% Channel inlet (x = 0) → v_in
for j = 1:M+2
    v_xP(j, 1) = v_in;
    v_yP(j, 1) = 0;
end

for i = 2:N+1
    for j = 2:M+1 

    %Blocking off method
    rr_re(j,i) = (x_P(i+1) - x_P(i)) / (((0.5 * (x_P(i+1) - x_P(i))) / (rho_ref / rho(j,i))) + ((0.5 * (x_P(i+1) - x_P(i))) / (rho_ref / rho(j,i+1))));
    rr_rw(j,i) = (x_P(i)-x_P(i-1)) / (((0.5 * (x_P(i)-x_P(i-1))) / (rho_ref / rho(j,i))) + ((0.5 * (x_P(i)-x_P(i-1))) / (rho_ref / rho(j,i-1))));
    rr_rn(j,i) = (y_P(j+1)-y_P(j)) / (((0.5 * (y_P(j+1)-y_P(j))) / (rho_ref / rho(j,i))) + ((0.5 * (y_P(j+1)-y_P(j))) / (rho_ref / rho(j+1,i))));
    rr_rs(j,i) = (y_P(j) - y_P(j-1)) / (((0.5 * (y_P(j) - y_P(j-1))) / (rho_ref / rho(j,i))) + ((0.5 * (y_P(j) - y_P(j-1))) / (rho_ref / rho(j-1,i))));
    
    % Velocity at faces         
    v_xn(j, i) = +(1/(rr_rn(j,i))) * (psi(j+1, i) - psi(j, i)) / (y_P(j+1)-y_P(j));
    v_xs(j, i) = +(1/(rr_rs(j,i))) * (psi(j, i) - psi(j-1, i)) / (y_P(j)-y_P(j-1));
    v_ye(j, i) = -(1/(rr_re(j,i))) * (psi(j, i+1) - psi(j, i)) / (x_P(i+1)-x_P(i));
    v_yw(j, i) = -(1/(rr_rw(j,i))) * (psi(j, i) - psi(j, i-1)) / (x_P(i)-x_P(i-1));
    
    % Velocity components at nodes
    v_xP(j, i) = (v_xn(j, i) + v_xs(j, i)) / 2;
    v_yP(j, i) = (v_ye(j, i) + v_yw(j, i)) / 2;

    end
end
    % Velocity magnitude
    v_P = sqrt(v_xP.^2 + v_yP.^2); 
end