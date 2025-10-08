function R = euler_angles_to_rotation_matrix(phi, theta, psi)
    % Compute rotation matrix from Euler angles (Z-Y-X convention)
    Rz = [cos(phi), -sin(phi), 0;
          sin(phi),  cos(phi), 0;
          0,         0,        1];

    Ry = [cos(theta), 0, sin(theta);
          0,          1, 0;
         -sin(theta), 0, cos(theta)];

    Rx = [1, 0,          0;
          0, cos(psi), -sin(psi);
          0, sin(psi),  cos(psi)];

    R = Rz * Ry * Rx;
end