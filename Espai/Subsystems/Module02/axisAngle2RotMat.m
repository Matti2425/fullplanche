function R = axisAngle2RotMat(axis,theta)
    % axis: [x; y; z] (3D vector)
    % theta: rotation angle (rad)
    % R: 3x3 rotation matrix

    % Normalizar eje
    axis = axis(:)/norm(axis);
    x = axis(1); y = axis(2); z = axis(3);

    c = cos(theta);
    s = sin(theta);

    R = [ x*x*(1-c) + c,     x*y*(1-c) - z*s, x*z*(1-c) + y*s;
          y*x*(1-c) + z*s,   y*y*(1-c) + c,   y*z*(1-c) - x*s;
          z*x*(1-c) - y*s,   z*y*(1-c) + x*s, c + z*z*(1-c)   ];
end