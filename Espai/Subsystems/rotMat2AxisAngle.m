function [axis, angle] = rotMat2AxisAngle(R)
    % rotMat2AxisAngle: Convierte matriz de rotación en eje-ángulo
    % R: matriz de rotación 3x3
    % axis: vector unitario [x; y; z]
    % angle: ángulo (rad)
    
    % Asegurar simetría numérica
    %R = (R + R')/2;
    
    % Calcular ángulo
    angle = acos((trace(R) - 1)/2);
    
    % Manejo de casos especiales
    if abs(angle) < 1e-12
        % Sin rotación
        axis = [1; 0; 0]; 
        angle = 0;
    elseif abs(angle - pi) < 1e-6
        % Rotación de 180 grados
        % El eje se extrae de la raíz de los elementos diagonales
        axis = sqrt((diag(R)+1)/2);
        % Ajustar signos
        if R(1,2) < 0, axis(2) = -axis(2); end
        if R(1,3) < 0, axis(3) = -axis(3); end
        axis = axis / norm(axis);
    else
        % Caso general
        axis = [R(3,2)-R(2,3);
                R(1,3)-R(3,1);
                R(2,1)-R(1,2)] / (2*sin(angle));
    end
end