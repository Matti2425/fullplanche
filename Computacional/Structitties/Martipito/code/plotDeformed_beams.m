function plotDeformed_beams(xn, Tn_st, u, scale, sigVM_st)
% Plots wing geometry (undeformed+deformed) and overlays stringers coloured by sigVM_st.
% Requires: plotDeformed.m (professor) in path, and wing.mat accessible for Tc in plotDeformed.

    % 1) Plot wing geometry using professor function (no stress)
    load('wing.mat','Tn_sk');              % uses skin as outer surface
    plotDeformed('wing', xn, Tn_sk, u, scale);  % <-- no sigVM
    hold on;

    % 2) Deformed nodal coordinates
    ux = u(1:6:end); uy = u(2:6:end); uz = u(3:6:end);
    x = xn(:,1) + scale*ux;
    y = xn(:,2) + scale*uy;
    z = xn(:,3) + scale*uz;

    % 3) Colour mapping (element-wise)
    c = sigVM_st(:);
    cmin = min(c); cmax = max(c);
    if abs(cmax-cmin) < 1e-14, cmax = cmin + 1; end
    cmap = jet(256);

    for e = 1:size(Tn_st,1)
        n1 = Tn_st(e,1); n2 = Tn_st(e,2);

        t = (c(e)-cmin)/(cmax-cmin);
        idx = 1 + round(t*(size(cmap,1)-1));
        idx = max(1, min(size(cmap,1), idx));

        plot3([x(n1) x(n2)], [y(n1) y(n2)], [z(n1) z(n2)], ...
              '-', 'Color', cmap(idx,:), 'LineWidth', 2.0);
    end

    colormap jet
    caxis([cmin cmax])  %Pa colorbar
    cb = colorbar;
    cb.Label.String = '\sigma_{VM} [Pa]';

    view(40,20);
    axis equal; axis vis3d;
end