function plotDeformed_beams(xn, Tn_st, u, scale, sigVM_st)
% plotDeformed_beams: deformed beam mesh coloured by element stress
% xn: [Nnodes x 3]
% Tn_st: [Nelem x 2]
% u: [Ndof x 1]
% scale: deformation scale
% sigVM_st: [Nelem x 1] (elementwise) or [] to colour by uz

    x0 = xn(:,1); y0 = xn(:,2); z0 = xn(:,3);

    ux = u(1:6:end); uy = u(2:6:end); uz = u(3:6:end);
    x = x0 + scale*ux;
    y = y0 + scale*uy;
    z = z0 + scale*uz;

    figure('Color','w'); hold on;

    % Undeformed discretised mesh (grey dashed)
    for e = 1:size(Tn_st,1)
        n1 = Tn_st(e,1); n2 = Tn_st(e,2);
        plot3([x0(n1) x0(n2)], [y0(n1) y0(n2)], [z0(n1) z0(n2)], ...
              ':', 'Color', 0.7*[1 1 1], 'LineWidth', 1.0);
    end

    % If stress not provided, colour by uz
    if nargin < 5 || isempty(sigVM_st)
        c = 0.5*(uz(Tn_st(:,1)) + uz(Tn_st(:,2)));   % element colour
        cb_label = 'u_z (arb.)';
    else
        c = sigVM_st(:);
        cb_label = '\sigma_{VM} [Pa]';
    end

    % Colormap mapping
    cmap = jet(256);
    cmin = min(c); 
    cmax = max(c);

    % Deformed mesh coloured
    for e = 1:size(Tn_st,1)
        n1 = Tn_st(e,1); n2 = Tn_st(e,2);
        t = (c(e)-cmin)/(cmax-cmin);
        idx = 1 + round(t*(size(cmap,1)-1));
        idx = max(1, min(size(cmap,1), idx));
        plot3([x(n1) x(n2)], [y(n1) y(n2)], [z(n1) z(n2)], ...
              '-', 'Color', cmap(idx,:), 'LineWidth', 2.0);
    end

    axis equal; axis vis3d; grid on;
    view(40,20);
    colormap jet
    caxis([cmin cmax]) %Pa
    cb = colorbar;
    cb.Label.String = '\sigma_{VM} [Pa]';
end
