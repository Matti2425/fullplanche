function eul = rotMat2Euler(R)
    % eul: [phi theta psi] (roll-pitch-yaw en rad)
    % Convención usada: ZYX (yaw-pitch-roll)
    
    % Extracción de ángulos según convención ZYX
    if abs(R(3,1)) < 1
        thetaY = -asin(R(3,1));
        thetaX = atan2(R(3,2), R(3,3));
        thetaZ = atan2(R(2,1), R(1,1));
    else
        % Singularidad (gimbal lock para ZYX)
        thetaY = -sign(R(3,1))*pi/2;
        thetaX = atan2(-R(1,2), -R(1,3));
        thetaZ = 0;
    end

    eul = [thetaX, thetaY, thetaZ]; % roll, pitch, yaw
end