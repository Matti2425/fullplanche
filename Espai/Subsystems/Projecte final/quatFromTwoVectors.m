function q = quatFromTwoVectors(u,v)
% Devuelve cuaternión scalar-first que rota u -> v
u = u(:) / norm(u);
v = v(:) / norm(v);
c = dot(u,v);
if c > 1 - 1e-12
    q = [1;0;0;0];
    return;
elseif c < -1 + 1e-12
    % antiparalelo: elegir un eje ortogonal arbitrario
    % Encuentra vector ortogonal a u
    if abs(u(1)) < abs(u(2))
        axis = [1;0;0];
    else
        axis = [0;1;0];
    end
    rot_axis = cross(u, axis);
    rot_axis = rot_axis / norm(rot_axis);
    q = [0; rot_axis]; % 180 deg rotation (scalar=0)
    return;
else
    w = cross(u,v);
    q = [1 + c; w];
    q = q / norm(q);
    return;
end
end
