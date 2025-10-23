function generate_Animation(r_POI_eci,r_sat_eci, R_E,t,n_steps)
%% --- Animation con coastlines rotando con w_E ---

figure('Color','w');
hold on; grid on; axis equal;

h_sat = plot3(r_sat_eci(1,1), r_sat_eci(2,1), r_sat_eci(3,1), 'ob', ...
              'MarkerFaceColor','b','MarkerSize',8);
h_poi = plot3(r_POI_eci(1,1), r_POI_eci(2,1), r_POI_eci(3,1), 'or', ...
              'MarkerFaceColor','r','MarkerSize',8);
h_vec = quiver3(0,0,0,0,0,0,'m','LineWidth',2,'MaxHeadSize',1);

title('Satellite Tracking Animation with rotating Earth');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
view(45,25)

[xs,ys,zs] = sphere(80);
surf(R_E*xs, R_E*ys, R_E*zs, ...
     'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.95);

%% --- Cargar coastlines ---
load coastlines
coast_lat = deg2rad(coastlat);
coast_lon = deg2rad(coastlon);
x_coast0 = R_E * cos(coast_lat) .* cos(coast_lon);
y_coast0 = R_E * cos(coast_lat) .* sin(coast_lon);
z_coast0 = R_E * sin(coast_lat);

% Dibujar la primera posición de las costas
h_coast = plot3(x_coast0, y_coast0, z_coast0, 'k', 'LineWidth', 1.1);

%% --- Parámetros ---
omega_E = 2*pi/86400;   % rad/s, rotación terrestre
dt = t(2) - t(1);

%% --- Animación ---
for k = 1:n_steps
    % --- Actualizar satélite y POI ---
    set(h_sat, 'XData', r_sat_eci(1,k), 'YData', r_sat_eci(2,k), 'ZData', r_sat_eci(3,k));
    set(h_poi, 'XData', r_POI_eci(1,k), 'YData', r_POI_eci(2,k), 'ZData', r_POI_eci(3,k));

    % --- Vector sat->POI ---
    vec = r_POI_eci(:,k) - r_sat_eci(:,k);
    set(h_vec, 'XData', r_sat_eci(1,k), 'YData', r_sat_eci(2,k), 'ZData', r_sat_eci(3,k), ...
               'UData', vec(1), 'VData', vec(2), 'WData', vec(3));

    % --- Rotar las coastlines con la Tierra ---
    theta_E = omega_E * t(k);
    Rz = [ cos(theta_E) -sin(theta_E) 0;
           sin(theta_E)  cos(theta_E) 0;
           0             0            1 ];

    r_coast_rot = Rz * [x_coast0'; y_coast0'; z_coast0'];
    set(h_coast, 'XData', r_coast_rot(1,:), 'YData', r_coast_rot(2,:), 'ZData', r_coast_rot(3,:));

    % --- Actualizar título ---
    title(sprintf('Satellite Tracking (step %d/%d) - Earth rotation: %.2f°', ...
          k, n_steps, rad2deg(theta_E)));

    drawnow;
end
end