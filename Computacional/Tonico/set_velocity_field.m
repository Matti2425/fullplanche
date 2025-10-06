function [u_P,v_P] = set_velocity_field(N, M, x_P, y_P)

% Computes the analytical velocity field for the Smith-Hutton problem
% at the cell centers of a 2D structured mesh.
%
%   u(x,y) =  2y(1 - x²)
%   v(x,y) = -2x(1 - y²)
%
% Inputs:
%   - N, M   : Number of divisions in x and y directions
%   - x_P, y_P : Coordinates of cell centers in x and y directions
%
% Outputs:
%   - u_P, v_P : 2D velocity field components at cell centers

u_P=zeros(M+2,N+2);
v_P=zeros(M+2,N+2);

for i=1:N+2
    for j=1:M+2
u_P(j,i)=2*y_P(j)*(1-x_P(i)^2);
v_P(j,i)=-2*x_P(i)*(1-y_P(j)^2);
    end
end

end

