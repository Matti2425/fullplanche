function q = quatMultiply(q1,q2)
% Producto de cuaterniones (scalar-first)
q0 = q1(1); v1 = q1(2:4);
r0 = q2(1); v2 = q2(2:4);
q = [ q0*r0 - dot(v1,v2);
      q0*v2 + r0*v1 + cross(v1,v2) ];
end