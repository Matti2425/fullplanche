function M = euler2RotMat(phi,theta,psi)
    % euler2DCM: Construye un DCM (ZYX) a partir de ángulos de Euler
    % phi   = roll  (deg) -> rotación sobre eje X
    % theta = pitch (deg) -> rotación sobre eje Y
    % psi   = yaw   (deg) -> rotación sobre eje Z
    % C     = DCM (3x3)

    phi = deg2rad(phi);
    theta = deg2rad(theta);
    psi = deg2rad(psi);

    % Rotación elemental sobre X (roll)
    Rx = [1     0          0;
          0 cos(phi) -sin(phi);
          0 sin(phi)  cos(phi)];

    % Rotación elemental sobre Y (pitch)
    Ry = [ cos(theta)  0  sin(theta);
           0           1  0;
          -sin(theta)  0  cos(theta)];

    % Rotación elemental sobre Z (yaw)
    Rz = [cos(psi) -sin(psi)  0;
          sin(psi)  cos(psi)  0;
          0         0         1];

    % Ensamblado final ZYX (estaba al revés)
    % M = Rz * Ry * Rx;
    M = Rx * Ry * Rz; 
end