% low_thrust_plane_r_polar.m
% Plantilla para optimización de transferencia de baja impulsión en plano
% Variables de estado: x = [r; theta; vr; h; m]
% Control: u = [T; alpha]  donde u_r = cos(alpha), u_theta = sin(alpha)
% Transcription directa con RK4 y fmincon

function low_thrust_plane_r_polar()
    clear; close all; clc;

    %%% --------------- Parámetros físicos ----------------
    mu = 3.986004418e14;      % [m^3/s^2] Tierra
    Isp = 300;                % [s]
    g0  = 9.80665;            % [m/s^2]
    Tmax = 0.3;               % [N] Empuje máximo (ejemplo muy pequeño)

    %%% --------------- Parámetros de problema --------------
    % Estados iniciales (ejemplo: órbita circular LEO)
    r0 = 7000e3;              % [m]
    theta0 = 0;               % [rad]
    vr0 = 0;                  % [m/s]
    h0 = sqrt(mu*r0);         % momento angular específico para circular
    m0 = 500;                 % [kg]

    % Objetivo: alcanzar r_f y h_f (órbita circular más alta)
    rf = 42164e3;             % [m] ejemplo: GEO (sólo como ejemplo)
    hf = sqrt(mu*rf);

    % Horizonte temporal y discretización
    tf = 200*24*3600;         % [s] tiempo total (ejemplo: 200 días)
    N = 120;                  % número de nodos (ajusta según necesidad)
    dt = tf/(N-1);

    % Regularización para suavizar controles (opcional)
    alpha_smooth_eps = 1e-4;

    %%% --------------- Construcción del vector de decisión --------
    % orden: [states(:,1); states(:,2); ...; states(:,N); controls(:,1); ...; controls(:,N)]
    nstate = 5;
    ncontrol = 2;
    nz = N*(nstate + ncontrol);

    % Bounds
    zlb = -inf(nz,1);
    zub =  inf(nz,1);

    % Para cada nodo: estados
    for k=1:N
        idx_s = state_index(k,N,nstate,ncontrol);
        % r:
        zlb(idx_s(1)) = 1e3;        % r mínima [m]
        zub(idx_s(1)) = 1e8;        % r máxima [m]
        % theta: libre
        zlb(idx_s(2)) = -inf; zub(idx_s(2)) = inf;
        % vr:
        zlb(idx_s(3)) = -2e4; zub(idx_s(3)) = 2e4;
        % h:
        zlb(idx_s(4)) = -1e9; zub(idx_s(4)) = 1e9;
        % m:
        zlb(idx_s(5)) = 1.0;  zub(idx_s(5)) = m0;
    end

    % Para cada nodo: controles [T; alpha]
    for k=1:N
        idx_u = control_index(k,N,nstate,ncontrol);
        % T:
        zlb(idx_u(1)) = 0;
        zub(idx_u(1)) = Tmax;
        % alpha:
        zlb(idx_u(2)) = -pi;
        zub(idx_u(2)) = pi;
    end

    % Initial guess z0: linea recta entre estados iniciales y finales
    z0 = zeros(nz,1);
    for k=1:N
        s = (k-1)/(N-1);
        idx_s = state_index(k,N,nstate,ncontrol);
        r_guess = r0 + s*(rf-r0);
        theta_guess = theta0 + s*2*pi;   % guess rotación (arbitrario)
        vr_guess = 0;
        h_guess = h0 + s*(hf-h0);
        m_guess = m0 - s*50;             % asumir consumo moderado
        z0(idx_s) = [r_guess; theta_guess; vr_guess; h_guess; m_guess];
    end
    for k=1:N
        idx_u = control_index(k,N,nstate,ncontrol);
        z0(idx_u) = [Tmax*0.5; 0];
    end

    %%% --------------- Llamada a fmincon -----------------
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
        'MaxIterations',2000,'MaxFunctionEvaluations',5e6,'TolFun',1e-8,'TolCon',1e-6);

    problem.objective = @(z) objective_total_thrust(z,N,dt,nstate,ncontrol,alpha_smooth_eps);
    problem.x0 = z0;
    problem.lb = zlb;
    problem.ub = zub;
    problem.nonlcon = @(z) dynamics_constraints(z,mu,Isp,g0,N,dt,nstate,ncontrol,r0,theta0,vr0,h0,m0,rf,hf);
    problem.options = options;
    problem.solver = 'fmincon';

    fprintf('Lanzando optimización con N=%d nodos, tf=%.1f días...\n',N,tf/86400);
    [z_opt,fval,exitflag,output] = fmincon(problem);

    %%% --------------- Extraer solución y mostrar ----------------
    [states,controls] = unpack_z(z_opt,N,nstate,ncontrol);
    tgrid = linspace(0,tf,N);

    figure; subplot(2,2,1); plot(tgrid/86400, states(1,:)/1e3); xlabel('t [d]'); ylabel('r [km]'); grid on;
    subplot(2,2,2); plot(tgrid/86400, states(3,:)); xlabel('t [d]'); ylabel('v_r [m/s]'); grid on;
    subplot(2,2,3); plot(tgrid/86400, states(4,:)); xlabel('t [d]'); ylabel('h [m^2/s]'); grid on;
    subplot(2,2,4); plot(tgrid/86400, states(5,:)); xlabel('t [d]'); ylabel('m [kg]'); grid on;

    figure; plot(states(1,:).*cos(states(2,:))/1e3, states(1,:).*sin(states(2,:))/1e3); axis equal;
    title('Trayectoria en el plano (km)'); xlabel('x [km]'); ylabel('y [km]'); grid on;

    fprintf('Optimización finalizada (exitflag=%d). fval=%.6f (J = integral T dt)\n', exitflag, fval);

end

%%% ----------------- Funciones auxiliares -----------------
function J = objective_total_thrust(z,N,dt,nstate,ncontrol,alpha_eps)
    % Minimizar integral de T (aproximado por suma Riemann)
    J = 0;
    for k=1:N
        idx_u = control_index(k,N,nstate,ncontrol);
        T_k = z(idx_u(1));
        J = J + T_k*dt;
    end
    % Suavizado en alpha (opcional): penaliza cambios grandes entre pasos
    if alpha_eps>0
        for k=1:N-1
            idx_u_k = control_index(k,N,nstate,ncontrol);
            idx_u_k1 = control_index(k+1,N,nstate,ncontrol);
            alpha_k = z(idx_u_k(2));
            alpha_k1 = z(idx_u_k1(2));
            J = J + alpha_eps*(alpha_k1 - alpha_k)^2;
        end
    end
end

function [c,ceq] = dynamics_constraints(z,mu,Isp,g0,N,dt,nstate,ncontrol,r0,theta0,vr0,h0,m0,rf,hf)
    % No hay desigualdades
    c = [];
    ceq = [];

    [states,controls] = unpack_z(z,N,nstate,ncontrol);

    % Dinámica por RK4 entre nodos
    for k=1:N-1
        xk = states(:,k);
        uk = controls(:,k);
        x_next = rk4_step(@(x,u) dynamics(x,u,mu,Isp,g0), xk, uk, dt);
        ceq = [ceq; states(:,k+1) - x_next]; %#ok<AGROW>
    end

    % Condiciones iniciales fijas
    x1 = states(:,1);
    ceq = [ceq;
        x1(1) - r0;
        x1(2) - theta0;
        x1(3) - vr0;
        x1(4) - h0;
        x1(5) - m0];

    % Condiciones finales: imponer rf y hf (no impongo theta final)
    xN = states(:,end);
    ceq = [ceq; xN(1) - rf; xN(4) - hf];

    % Masa positiva ya impuesta por bounds
end

function xdot = dynamics(x,u,mu,Isp,g0)
    % x = [r; theta; vr; h; m]
    r = x(1); theta = x(2); vr = x(3); h = x(4); m = x(5);
    T = u(1); alpha = u(2);

    ur = cos(alpha);
    utheta = sin(alpha);

    dr = vr;
    dtheta = h / r^2;
    dvr = h^2 / r^3 - mu / r^2 + (T/m)*ur;
    dh = r * (T/m) * utheta;
    dm = -T/(Isp*g0);

    xdot = [dr; dtheta; dvr; dh; dm];
end

function x_next = rk4_step(f, xk, uk, dt)
    k1 = f(xk, uk);
    k2 = f(xk + 0.5*dt*k1, uk);
    k3 = f(xk + 0.5*dt*k2, uk);
    k4 = f(xk + dt*k3, uk);
    x_next = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function idx = state_index(k,N,nstate,ncontrol)
    % devuelve indices en z para el estado en el nodo k
    base = (k-1)*(nstate + ncontrol);
    idx = (1: nstate) + base;
end

function idx = control_index(k,N,nstate,ncontrol)
    base = (k-1)*(nstate + ncontrol);
    idx = (nstate+1 : nstate + ncontrol) + base;
end

function [states,controls] = unpack_z(z,N,nstate,ncontrol)
    states = zeros(nstate,N);
    controls = zeros(ncontrol,N);
    for k=1:N
        idx_s = state_index(k,N,nstate,ncontrol);
        idx_u = control_index(k,N,nstate,ncontrol);
        states(:,k) = z(idx_s);
        controls(:,k) = z(idx_u);
    end
end
