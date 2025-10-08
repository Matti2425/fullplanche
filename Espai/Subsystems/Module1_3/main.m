clear; clc; close all;

% Definición de los subcuerpos

Bodies = []; % Vector de cuerpos vacío
Bodies = [Bodies, Body("Cuerpo central", 200, [0; 0; 0], diag([50, 60, 70]))];
Bodies = [Bodies, Body("Panel izquierdo", 20, [-2; 0; 0], diag([5, 30, 35]))];
Bodies = [Bodies, Body("Panel derecho", 20, [ 2; 0; 0], diag([5, 30, 35]))];


% Calculamos CG global

total_mass = sum([Bodies.mass]);
r_cg_total = [0; 0; 0];

for i = 1:length(Bodies)
    r_cg_total = r_cg_total + Bodies(i).mass * Bodies(i).r_cg;
end
r_cg_total = r_cg_total / total_mass;

fprintf('Centro de gravedad global [m]: (%.3f, %.3f, %.3f)\n', ...
        r_cg_total(1), r_cg_total(2), r_cg_total(3));


% Calculamos del tensor de inercia total respecto al CG global

I_total = zeros(3,3);
for i = 1:length(Bodies)
    I_total = I_total + Bodies(i).inertiaAtRef(r_cg_total);
end

disp('Tensor de inercia total respecto al CG global [kg·m²]:');
disp(I_total);


