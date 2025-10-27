function O = OmegaMatrix(omega)
% Matriz Ω(ω) usada en qdot = 0.5 * Omega * q  (scalar-first)
wx = omega(1); wy = omega(2); wz = omega(3);
O = [ 0   -wx  -wy  -wz;
      wx   0    wz  -wy;
      wy  -wz   0    wx;
      wz   wy  -wx   0 ];
end
