classdef Body
    
    properties
        name            % Nombre
        mass            % kg
        r_cg            % Vector posición del centro de gravedad en el ref. global [3x1]
        I_local         % Tensor de inercia respecto al centro de gravedad local [3x3]
    end

    methods
        % Constructor
        function obj = Body(name, mass, r_cg, I_local)
            if nargin > 0
                obj.name = name;
                obj.mass = mass;
                obj.r_cg = r_cg(:);       % Asegura que sea columna
                obj.I_local = I_local;
            end
        end
        
        % ---------- Método: tensor trasladado al sistema de referencia global ----------
        function I_global = inertiaAtRef(obj, r_ref)
            % Calcula el tensor de inercia respecto a un punto de referencia (r_ref)
            % usando el teorema de ejes paralelos.
            
            d = obj.r_cg - r_ref;  % Vector de desplazamiento
%             d_skew = [  0   -d(3)  d(2);
%                        d(3)   0   -d(1);
%                       -d(2)  d(1)   0  ];
%           se usa implícitamente
            I_global = obj.I_local + obj.mass * ( (d.'*d)*eye(3) - (d*d.') );
        end
    end
end