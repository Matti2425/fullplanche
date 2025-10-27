function qc = quatConj(q)
% Conjugado / inverso para cuaternión unitario (scalar-first)
qc = [ q(1); -q(2); -q(3); -q(4) ];
end