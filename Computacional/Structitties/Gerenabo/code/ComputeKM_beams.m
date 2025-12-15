function [K, M, f] = ComputeKM_beams(xn, N_nodes, N_elements, Tn, Tm, E_vector, G_vector, rho_vector, j_orientations, A_vector, Iy_vector, Iz_vector, J_vector, ky_vector, kz_vector, kt_vector, Fe, Be, Qe)
    
    % --- 1. Inicialización y Pre-computación ---
    Ndof = 6 * N_nodes;
    
    % Pre-asignación para el ensamblaje disperso (Sparse)
    % Cada elemento viga 3D tiene 12 DOFs, así que la matriz local es 12x12 = 144 entradas.
    max_entries = N_elements * 144;
    I_idx = zeros(max_entries, 1); % Índices de fila
    J_idx = zeros(max_entries, 1); % Índices de columna
    K_val = zeros(max_entries, 1); % Valores de Rigidez
    M_val = zeros(max_entries, 1); % Valores de Masa
    
    % Inicialización del vector de fuerzas
    f = zeros(Ndof, 1);
    
    % Constantes de integración de Gauss (2 puntos) fuera del bucle
    Xi = [-1/sqrt(3); 1/sqrt(3)];
    w = [1; 1];
    
    % Convertir cargas puntuales directamente (Vectorización)
    % Fe: [Valor, Nodo, DOF]
    if ~isempty(Fe)
        indices_f = 6*(Fe(:,2)-1) + Fe(:,3);
        f = f + accumarray(indices_f, Fe(:,1), [Ndof, 1]);
    end
    
    % Pre-procesamiento de cargas distribuidas en nodos (Qe y Be)
    % Acumulamos en matrices nodales para acceso rápido
    Q_nodes = zeros(N_nodes, 6);
    if ~isempty(Qe)
        for r = 1:size(Qe, 1)
            Q_nodes(Qe(r,2), Qe(r,3)) = Q_nodes(Qe(r,2), Qe(r,3)) + Qe(r,1);
        end
    end
    
    B_nodes = zeros(N_nodes, 6);
    if ~isempty(Be)
        for s = 1:size(Be, 1)
            B_nodes(Be(s,2), Be(s,3)) = B_nodes(Be(s,2), Be(s,3)) + Be(s,1);
        end
    end

    % --- 2. Bucle sobre Elementos ---
    idx_count = 0; % Contador para llenar los vectores sparse
    
    for e = 1:N_elements
        % Nodos y Material del elemento actual
        n1 = Tn(e, 1);
        n2 = Tn(e, 2);
        mat_idx = Tm(e, 1);
        
        % Geometría
        p1 = xn(n1, :)';
        p2 = xn(n2, :)';
        vec_L = p2 - p1;
        Le = norm(vec_L);
        
        % Sistema de Coordenadas Local (Matriz de Rotación)
        u = vec_L / Le;           % Eje x local
        v_temp = j_orientations(mat_idx, :)';
        w_vec = cross(u, v_temp); % Eje z local
        w_vec = w_vec / norm(w_vec); % Normalizar por seguridad
        v_vec = cross(w_vec, u);  % Eje y local recalculado ortogonalmente
        
        % Matriz de rotación base (3x3)
        R0 = [u, v_vec, w_vec]'; 
        % Matriz de rotación total (12x12) usando bloque diagonal
        Re = blkdiag(R0, R0, R0, R0);
        
        % --- Matrices de Rigidez (Stiffness) ---
        % Propiedades
        E = E_vector(mat_idx); 
        G = G_vector(mat_idx);
        A = A_vector(mat_idx);
        Iy = Iy_vector(mat_idx); Iz = Iz_vector(mat_idx);
        Jt = J_vector(mat_idx);
        ky = ky_vector(mat_idx); kz = kz_vector(mat_idx);
        
        % Funciones de forma derivadas (constantes)
        invL = 1/Le;
        Nx_diff = [-invL, invL]; % dN/dx
        
        % Axial (Ka)
        Ba = zeros(1, 12);
        Ba([1, 7]) = Nx_diff;
        Ka = (Le * E * A) * (Ba' * Ba);
        
        % Bending (Kb)
        Bb = zeros(2, 12);
        Bb(1, [5, 11]) = Nx_diff; % Derivada para rotación y
        Bb(2, [6, 12]) = Nx_diff; % Derivada para rotación z
        Cb = diag([E*Iy, E*Iz]);
        Kb = Le * (Bb' * Cb * Bb);
        
        % Shear (Ks)
        N_Gauss = 0.5; % Valor medio para integración reducida en cortante
        Bs = zeros(2, 12);
        % Filas correspondientes a cortante en y y z
        Bs(1, [2, 6, 8, 12]) = [Nx_diff(1), -N_Gauss, Nx_diff(2), -N_Gauss]; 
        Bs(2, [3, 5, 9, 11]) = [Nx_diff(1),  N_Gauss, Nx_diff(2),  N_Gauss];
        Cs = G * A * diag([ky, kz]);
        Ks = Le * (Bs' * Cs * Bs);
        
        % Torsion (Kt)
        Bt = zeros(1, 12);
        Bt([4, 10]) = Nx_diff;
        Kt = (Le * G * Jt * kt_vector(mat_idx)) * (Bt' * Bt);
        
        % Rigidez Local Total (en coordenadas locales)
        Ke_local = Ka + Kb + Ks + Kt;
        
        % Rotar a coordenadas globales: K_global = R' * K_local * R
        % Nota: En tu código original hacías R'*B'*C*B*R. 
        % Matemáticamente (R'*B'*C*B*R) es lo mismo que R'*(K_local)*R
        Ke_global = Re' * Ke_local * Re;
        
        % --- Matriz de Masa (Mass) y Fuerzas Distribuidas ---
        Me_local = zeros(12, 12);
        fe_local_dist = zeros(12, 1);
        
        rho = rho_vector(mat_idx);
        Rho_mat = diag([rho*A, rho*A, rho*A, rho*Jt, rho*Iy, rho*Iz]);
        
        % Cargas en los nodos del elemento actual
        q_vec = [Q_nodes(n1, :), Q_nodes(n2, :)]';
        b_vec = [B_nodes(n1, :), B_nodes(n2, :)]';
        
        % Integración de Gauss
        for k = 1:2
            Ni = (1 + Xi(k)*[-1, 1]) / 2; % [N1, N2]
            
            % Matriz de funciones de forma N (6x12)
            % Se puede construir usando producto de Kronecker para eficiencia
            % N_mat = [N1*I, N2*I] donde I es identidad 6x6
            N_mat = kron(Ni, eye(6)); 
            
            % Masa
            detJ = Le / 2; % Determinante del Jacobiano
            Me_local = Me_local + (w(k) * detJ * (N_mat' * Rho_mat * N_mat));
            
            % Fuerzas Distribuidas (Integración)
            % fe = Integral(N' * N * q) + Masa * b
            % Nota: Tu código original multiplica por R al principio y al final.
            % Asumiré que q y b están en globales y deben rotarse o proyectarse.
            % Siguiendo lógica original:
            term_force = N_mat' * N_mat * Re * q_vec; 
            fe_local_dist = fe_local_dist + (w(k) * detJ * Re' * term_force);
        end
        
        % Rotar matriz de masa
        Me_global = Re' * Me_local * Re;
        
        % Añadir contribución de fuerzas de cuerpo (Me * b_vec)
        fe_total = fe_local_dist + Me_global * b_vec;
        
        % --- Preparar Índices para Ensamblaje ---
        dof_indices = [ (6*(n1-1)+1 : 6*n1), (6*(n2-1)+1 : 6*n2) ];
        
        % Generar rejilla de índices (meshgrid es rápido)
        [C_idx, R_idx] = meshgrid(dof_indices, dof_indices);
        
        % Llenar vectores de almacenamiento
        start_idx = idx_count + 1;
        end_idx = idx_count + 144;
        
        I_idx(start_idx:end_idx) = R_idx(:);
        J_idx(start_idx:end_idx) = C_idx(:);
        K_val(start_idx:end_idx) = Ke_global(:);
        M_val(start_idx:end_idx) = Me_global(:);
        
        idx_count = end_idx;
        
        % Ensamblar vector de fuerza directamente (es denso y rápido)
        f(dof_indices) = f(dof_indices) + fe_total;
    end
    
    % --- 3. Ensamblaje Final de Matrices ---
    % Sumar duplicados automáticamente con sparse
    K = sparse(I_idx, J_idx, K_val, Ndof, Ndof);
    M = sparse(I_idx, J_idx, M_val, Ndof, Ndof);
    
end