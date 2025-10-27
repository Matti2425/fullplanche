function generate_Animation_GPT(r_POI_eci, r_sat_eci, R_E, t, n_steps,w_E, w_sat)
% generate_Animation - Animación de la órbita, el POI y los ejes BRF del satélite
%
% Entradas:
%   r_POI_eci   : 3xN posiciones del POI en ECI
%   r_sat_eci   : 3xN posiciones del satélite en ECI
%   R_E         : radio terrestre (m)
%   t           : vector temporal (s)
%   n_steps     : pasos de animación
%   w_sat       : velocidad angular [3x1] o [3xN] del satélite (rad/s)

%% --- Configuración de la figura ---
figure('Color','w');
hold on; grid on; axis equal;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(45,25);
title('Satellite tracking with Earth and rotating BRF');

%% --- Dibujar la Tierra con las costas ---
% Superficie esférica
[xE, yE, zE] = sphere(100);
surf(R_E*xE, R_E*yE, R_E*zE, 'FaceColor', [0.6 0.8 1], ...
     'EdgeColor','none', 'FaceAlpha',0.7);

% Cargar y trazar las líneas de costa
load coastlines
% Convertir lat/lon (grados) a radianes y proyectarlas a ECI (Z hacia el polo norte)
coast_lat = deg2rad(coastlat);
coast_lon = deg2rad(coastlon);
x_coast0 = R_E * cos(coast_lat) .* cos(coast_lon);
y_coast0 = R_E * cos(coast_lat) .* sin(coast_lon);
z_coast0 = R_E * sin(coast_lat);
h_coast = plot3(x_coast0, y_coast0, z_coast0, 'k', 'LineWidth', 1.1);

%% --- Elementos dinámicos ---
h_sat = plot3(r_sat_eci(1,1), r_sat_eci(2,1), r_sat_eci(3,1), 'ob', ...
              'MarkerFaceColor','b','MarkerSize',8);
h_poi = plot3(r_POI_eci(1,1), r_POI_eci(2,1), r_POI_eci(3,1), 'or', ...
              'MarkerFaceColor','r','MarkerSize',8);

% Vector sat->POI
v = r_POI_eci(:,1) - r_sat_eci(:,1);
h_vec = quiver3(r_sat_eci(1,1), r_sat_eci(2,1), r_sat_eci(3,1), ...
                v(1), v(2), v(3), 0, 'm', 'LineWidth',2,'MaxHeadSize',1);

% Ejes BRF (body)
h_bx = plot3([0 0],[0 0],[0 0],'-','LineWidth',2,'Color',[0.85 0.33 0.10]); % rojo
h_by = plot3([0 0],[0 0],[0 0],'-','LineWidth',2,'Color',[0.47 0.67 0.19]); % verde
h_bz = plot3([0 0],[0 0],[0 0],'-','LineWidth',2,'Color',[0 0.4470 0.7410]); % azul

%% --- Inicialización de rotación ---
R_bi = eye(3); % orientación inicial identidad
N = length(t);
dt = mean(diff(t));
if isvector(w_sat) && numel(w_sat)==3
    w_all = repmat(w_sat(:),1,N);
elseif size(w_sat,2)==N
    w_all = w_sat;
else
    w_all = zeros(3,N);
end

axis_len = 2e6; % longitud de ejes BRF (ajusta si se ve grande/pequeño)

%% --- Animación principal ---
for k = 1:min(N, n_steps)
    rs = r_sat_eci(:,k);
    rp = r_POI_eci(:,k);

    % Actualizar satélite y POI
    set(h_sat, 'XData', rs(1), 'YData', rs(2), 'ZData', rs(3));
    set(h_poi, 'XData', rp(1), 'YData', rp(2), 'ZData', rp(3));

    % Vector sat->POI
    v = rp - rs;
    delete(h_vec);
    h_vec = quiver3(rs(1), rs(2), rs(3), v(1), v(2), v(3), 0, ...
                    'm','LineWidth',2,'MaxHeadSize',1);

     % --- Rotar las coastlines con la Tierra ---
    theta_E = w_E * t(k);
    Rz = [ cos(theta_E) -sin(theta_E) 0;
           sin(theta_E)  cos(theta_E) 0;
           0             0            1 ];

    r_coast_rot = Rz * [x_coast0'; y_coast0'; z_coast0'];
    set(h_coast, 'XData', r_coast_rot(1,:), 'YData', r_coast_rot(2,:), 'ZData', r_coast_rot(3,:));

    % --- Actualizar BRF ---
    omega = w_all(:,k);
    R_bi = R_bi * expm(skew(omega)*dt); % rotación incremental

    bx = R_bi(:,1); by = R_bi(:,2); bz = R_bi(:,3);

    set(h_bx, 'XData',[rs(1), rs(1)+bx(1)*axis_len], ...
              'YData',[rs(2), rs(2)+bx(2)*axis_len], ...
              'ZData',[rs(3), rs(3)+bx(3)*axis_len]);
    set(h_by, 'XData',[rs(1), rs(1)+by(1)*axis_len], ...
              'YData',[rs(2), rs(2)+by(2)*axis_len], ...
              'ZData',[rs(3), rs(3)+by(3)*axis_len]);
    set(h_bz, 'XData',[rs(1), rs(1)+bz(1)*axis_len], ...
              'YData',[rs(2), rs(2)+bz(2)*axis_len], ...
              'ZData',[rs(3), rs(3)+bz(3)*axis_len]);

    title(sprintf('Satellite tracking – t = %.1f s', t(k)));
    drawnow limitrate
end

end
