function spherical_plot(r, theta, phi)
    % Función para graficar puntos o una superficie en coordenadas esféricas
    % r: Radio de la esfera o vector de radios
    % theta: Ángulo azimutal en radianes (vector)
    % phi: Ángulo polar en radianes (vector)
    
    % Comprobamos si r es un valor escalar o un vector
    if numel(r) == 1
        % Si r es un único valor (radio constante), lo replicamos para todos los puntos
        r = r * ones(size(theta));
    end
    
    % Convertir coordenadas esféricas a coordenadas cartesianas
    X = r .* sin(phi) .* cos(theta);
    Y = r .* sin(phi) .* sin(theta);
    Z = r .* cos(phi);
    
    % Crear el gráfico de superficie o puntos
    figure;
    plot3(X, Y, Z); % Usamos surf para dibujar una superficie
    % Si prefieres un gráfico de puntos, usa plot3(X, Y, Z) en su lugar
    
    % Mejora de la visualización
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Gráfico Esférico');
    
    % Agregar una vista 3D adecuada
    view(3);
end