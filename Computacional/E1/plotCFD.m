function plotCFD(titleStr, F, N, M, dx, dy)
% plotCFD generates a graphical representation of the scalar field F using a color map.
%
% Inputs:
%   titleStr - Title of the graph (string).
%   F        - Matrix representing the scalar field to be plotted.
%   N        - Number of discretization points in the x-direction.
%   M        - Number of discretization points in the y-direction.
%   dx       - Grid spacing in the x-direction.
%   dy       - Grid spacing in the y-direction.
%
% Outputs:
%   A graphical representation of F with color coding and a color bar.
%--------------------------------------------------------------------------

% Plot the field F
imagesc((0.5:N-0.5) * dx, (0.5:M-0.5) * dy, F(2:end-1, 2:end-1)); 
axis xy equal tight;
colorbar;
shading interp; 
colormap jet;

% Set axis labels
xlabel('x [m]'); 
ylabel('y [m]'); 
grid off;

% Set the title correctly
if isempty(titleStr)
    title('Function F', 'FontWeight', 'bold', 'FontSize', 12);
else
    title(['Function ', titleStr], 'FontWeight', 'bold', 'FontSize', 12);
end

end