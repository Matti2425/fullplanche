function plotResults(surf_x, surf_y, phi_def, nod_x, nod_y, Geometry, Fluid, diag, case_study)
% PLOTRESULTS Generates visualization plots for the convection-diffusion solution
%
% This function creates surface plots of the phi field for different gamma values
% and optionally compares with analytical solutions for the Smith-Hutton case.
%
% INPUTS:
%   surf_x, surf_y - Meshgrid coordinates
%   psi_def        - Solution array (N+2, M+2, num_gamma)
%   nod_x, nod_y   - Node coordinates
%   Geometry       - Structure with H
%   Fluid          - Structure with gamma
%   diag           - Bool: true for diagonal, false for Smith-Hutton
%   case_study     - case_study selected
%--------------------------------------------------------------------------

    gamma = Fluid.gamma;
    H = Geometry.H;
    
    if ~diag
        % Smith-Hutton case: Compare with analytical solution
        r = [0.000, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000]';
        p = [1.989, 1.402, 1.146, 0.946, 0.775, 0.621, 0.480, 0.349, 0.227, 0.111, 0.000;
             2.00,  1.999, 1.9997, 1.985, 1.841, 0.951, 0.154, 0.001, 0.000, 0.000, 0.000;
             2.000, 2.000, 2.000, 1.999, 1.964, 1.000, 0.036, 0.001, 0.000, 0.000, 0.000]';
        
        phi_n = zeros(length(r), size(p, 2));
        
        % Interpolate numerical solution at analytical points
        for i = 1:size(p, 2)
            assert(r(1) >= (nod_x(1) - H), 'You must have min(x) <= min(xx).');
            assert(r(end) <= (nod_x(end) - H), 'You must have max(xx) <= max(x).');
            
            int_phi = interp1(nod_x - H, phi_def(:, 1, i), r, 'linear', 'extrap');
            phi_n(:, i) = int_phi(:);
        end
        
        % Verification plot
        figure('Name', 'Code Verification', 'NumberTitle', 'off');
        hold on;
        for i = 1:length(gamma)
            plot(r, phi_n(:, i), 'LineWidth', 2);
            plot(r, p(:, i), 'o', 'MarkerSize', 5);
        end
        grid on;
        xlabel('x [m]');
        ylabel('\phi');
        ylim([-0.1, 2.1]);
        s = size(surf_x) - 2;
        mesh_size = s(1) + "x" + s(2);
        name = "Code Verification (Analytic vs Numeric) " + mesh_size;
        title(name, 'FontSize', 14);
        legend('\rho/\Gamma = 10', 'Analytic \rho/\Gamma = 10', ...
               '\rho/\Gamma = 10^3', 'Analytic \rho/\Gamma = 10^3', ...
               '\rho/\Gamma = 10^6', 'Analytic \rho/\Gamma = 10^6');
        hold off;
        dir = "figures/errors/Code Verification_" + mesh_size + ".png";
        exportgraphics(gca, dir, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'white');
    end
    
    % Distribution maps for each gamma value
    for i = 1:length(gamma)
        figure('Name', sprintf('Distribution Map - rho/Gamma = %.0e', Fluid.GAMMA(i)), 'NumberTitle', 'off');
        surf(surf_x, surf_y, phi_def(:, :, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(0, 90);
        c = colorbar;
        c.Label.String = '\Psi';
        colormap jet;
        title(sprintf('Distribution Map \\rho/\\Gamma = %.0e', Fluid.GAMMA(i)));
        % Texto tipo nota arriba a la derecha con fondo blanco
        annotation('textbox', [0.69, 0.75, 0.25, 0.1], ...
            'String', "Mesh Size: " + mesh_size, ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 12);


        axis equal;
        xlabel('x [m]');
        ylabel('y [m]');
        grid off;
        
        if ~diag
            xlim([-1, 1]);
            ylim([0, 1]);
        else
            xlim([-10, 0]);
            ylim([0, 10]);
        end
        dir = "figures/" + case_study + "/Distribution_map_" + Fluid.GAMMA(i) + "_" + mesh_size + ".png";
        % Obtener la carpeta de la ruta
        folder = fileparts(dir);
        
        % Crear la carpeta si no existe
        if ~exist(folder, "dir")
            mkdir(folder);
        end
        exportgraphics(gca, dir, 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'white');
    end

end