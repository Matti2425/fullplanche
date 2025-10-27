function S = skew(v)
% skew - devuelve la matriz cruz (o matriz skew-symmetric) de un vector 3x1
%
%   S * a  ==  cross(v, a)
%
%   Entrada:
%       v : vector 3x1
%   Salida:
%       S : matriz 3x3 tal que S*a = v x a

S = [   0   -v(3)   v(2);
       v(3)    0   -v(1);
      -v(2)  v(1)    0  ];
end
