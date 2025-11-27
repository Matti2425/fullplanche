function plotOrbit(a, style, labelStr, varargin)
    theta = linspace(0,2*pi,400);
    x = a*cos(theta); y = a*sin(theta);
    if nargin > 3
        % transfer ellipse
        r1 = varargin{1}; r2 = varargin{2};
        a_t = a; e_t = abs(r2 - r1)/(r2 + r1);
        r = a_t*(1 - e_t^2)./(1 + e_t*cos(theta));
        x = r.*cos(theta); y = r.*sin(theta);
    end
    plot(x, y, style, 'DisplayName', labelStr, 'LineWidth',1.2);
end