function plotDeformed(file, xn, Tn, u, scale, sigVM)
    % Load data
    switch file
        case 'shell'
            load('shell.mat', 'Tc');
        case 'wing'
            load('wing.mat', 'Tc');
    end

    % Precompute
    x0 = xn(:, 1);
    y0 = xn(:, 2);
    z0 = xn(:, 3);
    ux = u(1:6:end);
    uy = u(2:6:end);
    uz = u(3:6:end);

    % Deformed coordinates
    x = x0 + scale * ux;
    y = y0 + scale * uy;
    z = z0 + scale * uz;

    % Initialize figure
    figure('Color', 'w');
    hold on;

    % Plot undeformed contour (wireframe)
    for e = 1:size(Tn, 1)
        % Extract the node indices for each face from Tn
        nodes = Tn(e, :);
        
        % Plot the wireframe for undeformed shape
        patch(x0(nodes), y0(nodes), z0(nodes), ones(1, 4), 'FaceColor', 'none', 'EdgeColor', 0.5 * [1, 1, 1]);
    end

    % Scale the displacement values for better color contrast
    uz_scaled = scale * uz; % Apply scaling to displacement to increase visual effect

    % Plot deformed structure (colored by Z displacement)
    for e = 1:size(Tn, 1)
        nodes = Tn(e, :);
        
        % Get the Z-displacement values for each node of the face
        cdata_face = uz_scaled(nodes); % Use scaled displacement for coloring
        
        % Plot the deformed shape, using the Z-displacement for color
        p = patch(x(nodes), y(nodes), z(nodes), cdata_face, ...
                  'FaceColor', 'interp', 'EdgeColor', 'k', 'FaceAlpha', 0.9);
    end

    % Option to plot Von Mises stress (if available)
    if exist('sigVM', 'var')
        SigVM = InterpFunction(sigVM, Tn);  % Interpolation of the Von Mises stress
        for e = 1:size(Tn, 1)
            nodes = Tn(e, :);
            patch(x(nodes), y(nodes), z(nodes), SigVM(e), ...
                  'FaceColor', 'interp', 'EdgeColor', 'none');
        end
        % Color axis for stress values
        clim([0, max([SigVM(:); 1])]); 
        cb = colorbar;
        set(cb, 'Ticks', [0, max([SigVM(:); 1])]);
    end

    % Additional options
    title(sprintf('Deformation Scale = %.1f', scale));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; axis vis3d;
    grid on;

    % Color bar for Z-displacement
    colorbar; colormap jet;
    caxis([min(uz_scaled), max(uz_scaled)]); % Apply caxis for better color mapping
    view(3);

end

function sign = InterpFunction(sige, Tn)
    a = [-1, 1, 1, -1];
    b = [-1, -1, 1, 1];
    sign = zeros(size(Tn, 1), 1);
    for e = 1:size(Tn, 1)
        N4 = zeros(4, 4);
        for k = 1:4
            for i = 1:4
                N4(k, i) = (1 + a(i) * a(k) / sqrt(3)) * (1 + b(i) * b(k) / sqrt(3)) / 4;
            end
        end
        sign(:, e) = N4 \ (sige(e, :)');
    end
end
