%% Simulación de la rotación libre de un cuerpo rígido
% Ecuaciones de Euler - efecto Dzhanibekov

clear; clc; close all

% Momentos principales de inercia
I1 = 1;      % eje menor
I2 = 2;      % eje intermedio
I3 = 3;      % eje mayor

% Condiciones iniciales
% Rotación casi alrededor del eje intermedio (inestable)
omega0 = [0.01; 1; 0.02];   % [ω1, ω2, ω3] iniciales

% Intervalo temporal
tspan = [0 50];

% Integración numérica
[t, w] = ode45(@(t,omega) euler_equations(t,omega,I1,I2,I3), tspan, omega0);

% Energía y momento angular
T = 0.5 * (I1*w(:,1).^2 + I2*w(:,2).^2 + I3*w(:,3).^2);
L2 = (I1*w(:,1)).^2 + (I2*w(:,2)).^2 + (I3*w(:,3)).^2;

% Graficas
figure
plot3(w(:,1), w(:,2), w(:,3), 'b', 'LineWidth', 1.5)
xlabel('\omega_1'); ylabel('\omega_2'); zlabel('\omega_3');
title('Trayectoria del vector de velocidad angular en el sistema del cuerpo')
grid on; axis equal

figure
plot(t, w)
xlabel('Tiempo'); ylabel('\omega_i (rad/s)')
legend('\omega_1','\omega_2','\omega_3')
title('Evolución temporal de las componentes de velocidad angular')

figure
plot(t, T, t, L2)
xlabel('Tiempo'); ylabel('Magnitudes conservadas')
legend('Energía cinética','|L|^2')
title('Conservación de energía y momento angular')
grid on

%% --- Función ecuaciones de Euler ---
function domega = euler_equations(~, omega, I1, I2, I3)
    domega = zeros(3,1);
    domega(1) = ((I2 - I3)/I1) * omega(2) * omega(3);
    domega(2) = ((I3 - I1)/I2) * omega(3) * omega(1);
    domega(3) = ((I1 - I2)/I3) * omega(1) * omega(2);
end
