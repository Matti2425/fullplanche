function q = axisAngle2Quat(axis,theta)
    % axis: [x; y; z] (3D vector)
    % theta: rotation angle (rad)
    % q: quaternion [w x y z] (w=scalar part)

    axis = axis(:)/norm(axis);
    s = sin(theta/2);

    q = [cos(theta/2);
         axis(1)*s;
         axis(2)*s;
         axis(3)*s];
end