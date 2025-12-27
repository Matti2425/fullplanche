% m2_hall_thruster.m
% Hall Effect Thruster - Modelo 1D no dimensional
% Basado estrictamente en M6_Lecture1
% ---------------------------------------------------------

clear; close all; clc

%% --- CONSTANTES FÍSICAS ---
qe = 1.6022e-19;
me = 9.1034e-31;
mi = 2.1802e-25;
kB = 1.380649e-23;

%% --- DATOS DEL THRUSTER ---
H  = 15e-3;
rm = 50e-3;
A  = 2*pi*rm*H;

cn = 300;
Gamma_n0 = 6e22;
B0 = 0.02;
L  = 10e-3;

Vi = 12.1;
Ei = qe*Vi;
Ei_p = 2.5*Ei;

sigma0 = 3.6e-20;
sigma_en = 27e-20;
alphaB = 1/16;

%% --- NORMALIZACIÓN ---
Tstar = Ei;
vstar = sqrt(Tstar/mi);
lstar = 0.01;
nstar = 1/(sigma0*lstar);
Gstar = nstar*vstar;
nustar = vstar/lstar;

%% --- PARÁMETROS DEL ESTUDIO ---
kTeB_over_Ei = [0.08 0.10 0.12 0.14 0.16];

%% --- CONDICIÓN DE ENTRADA ---
Gamma_d_hat = -0.01;

%% --- MALLA AXIAL ---
xspan = [0 1];

%% --- ALMACENAMIENTO ---
SOL = struct();

%% === BUCLE PRINCIPAL ===
for k = 1:length(kTeB_over_Ei)

    TeB_hat = kTeB_over_Ei(k);

    % Condiciones iniciales (x = 0)
    vi0   = 0.05;
    ne0   = 0.01;
    Te0   = TeB_hat;
    phi0  = 0;

    y0 = [vi0; ne0; Te0; phi0];

    % Integración ODE
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    sol = ode15s(@(x,y) HET_ODE(x,y,Gamma_d_hat,...
                cn,alphaB,B0,...
                sigma0,sigma_en,...
                Ei_p,me,mi,kB), xspan, y0, opts);

    x = linspace(0,1,400);
    Y = deval(sol,x);

    vi  = Y(1,:);
    ne  = Y(2,:);
    Te  = Y(3,:);
    phi = Y(4,:);

    M = vi ./ sqrt(Te);

    % --- EFICIENCIAS ---
    ci = sqrt(2*qe*abs(phi(end))*Tstar/mi);
    nn = (Gamma_n0 - ne.*vi) ./ cn;

    eta_u = (ne(end)*vi(end)) / (nn(1)*cn);
    eta_a = 1 - abs(Gamma_d_hat)/(ne(end)*vi(end));

    SOL(k).x = x;
    SOL(k).vi = vi;
    SOL(k).ne = ne;
    SOL(k).Te = Te;
    SOL(k).phi = phi;
    SOL(k).M = M;
    SOL(k).eta_u = eta_u;
    SOL(k).eta_a = eta_a;
end

%% === GRÁFICAS (PREGUNTA 1) ===
figure
for k=1:length(SOL)
    plot(SOL(k).x, SOL(k).ne,'LineWidth',1.5); hold on
end
xlabel('x/L'), ylabel('n_e')
legend('0.08','0.10','0.12','0.14','0.16')
title('Electron density')

figure
for k=1:length(SOL)
    plot(SOL(k).x, SOL(k).M,'LineWidth',1.5); hold on
end
xlabel('x/L'), ylabel('Mach')
title('Mach number distribution')

%% === EFICIENCIAS (PREGUNTA 2) ===
eta_u = [SOL.eta_u];
eta_a = [SOL.eta_a];

figure
plot(kTeB_over_Ei,eta_u,'o-',kTeB_over_Ei,eta_a,'s-','LineWidth',1.5)
xlabel('kT_{eB}/E_i'), ylabel('\eta')
legend('\eta_u','\eta_a')
title('Utilization and back-streaming efficiency')

%% === CASO B(x) NO UNIFORME (PREGUNTA 3) ===
Bx = @(x) 0.04*exp(-((x-0.5).^2)/(2*(0.1)^2));

figure
plot(linspace(0,1,400),Bx(linspace(0,1,400)),'LineWidth',2)
xlabel('x/L'), ylabel('B(x) [T]')
title('Non-uniform magnetic field')

%% Funciones auxiliares
function dydx = HET_ODE(x, y, Gamma_d_hat, ...
                        cn, alphaB, B, ...
                        sigma0, sigma_en, ...
                        Ei_p, me, mi, kB)
% Sistema ODE 1D no dimensional para Hall Effect Thruster
% Basado en M6_Lecture1
%
% Variables:
%   y(1) = v_i   velocidad iónica adimensional
%   y(2) = n_e   densidad electrónica adimensional
%   y(3) = T_e   temperatura electrónica adimensional
%   y(4) = phi   potencial adimensional

vi  = y(1);
ne  = max(y(2),1e-8);   % evitar valores no físicos
Te  = max(y(3),1e-8);
phi = y(4);

nn = max((1 - ne.*vi)/cn, 1e-8);   % conservación de flujo neutral

% Freq ionización
Rion = sigma0 * sqrt(Te);          % ley simplificada de ionización
nu_ion = nn .* Rion;

% Fuente
Sion = nu_ion .* ne;

% Bohm
omega_ce = abs(B)*1.6022e-19/me;
nu_e = alphaB * omega_ce;
mu_e = 1 ./ nu_e;

% EDOs
% 1) Continuidad iónica
dnedx = (Sion - ne .* gradient_safe(vi)) ./ max(vi,1e-6);

% 2) Momento iónico
dvidx = -gradient_safe(phi) ./ max(vi,1e-6);

% 3) Energía electrónica
dTedx = (2/3) * ( ...
         nu_ion .* Ei_p ...
         - sigma_en .* nn .* Te.^(3/2) );

% 4) Corriente electrónica
dphidx = -(Gamma_d_hat) ./ (mu_e .* ne);

% Vector salida
dydx = [dvidx;
        dnedx;
        dTedx;
        dphidx];
end

% Gradiente
function g = gradient_safe(f)
% Evita singularidades en derivadas
g = f;
g(abs(g) < 1e-6) = sign(g(abs(g) < 1e-6))*1e-6;
end

