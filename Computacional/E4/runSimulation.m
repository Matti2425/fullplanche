function results = runSimulation(params)
num_cases = numel(params.cases);
results = repmat(struct("type", [], "N", [], "Ck", [], "k", [], "E_k", [], "time", [], "energy", []), num_cases, 1);

for idx = 1:num_cases
    case_info = params.cases(idx);
    N = case_info.N;
    case_type = case_info.type;
    
    fprintf("  Case %d/%d: %s N=%d", idx, num_cases, case_type, N);
    if strcmp(case_type, 'LES')
        fprintf(" Ck=%.4f", case_info.Ck);
    end
    fprintf("\n");
    
    % Initialize: use positive modes only 
    k = 1:N;  % Mode numbers
    u = 1 ./ k;  % Initial velocities
    time_step = params.C1 * params.Re / (N^2);
    
    % Iterative solver
    for iter = 1:params.maxIt
        % Preallocate
        conv = zeros(1, N);
        diff = zeros(1, N);
        
        % Compute convective and diffusive terms for each mode j
        for j = 1:N
            % Convective term 
            for p = (-N+j):N
                q = j - p;
                if p ~= 0 && q ~= 0
                    % Handle positive and negative modes (matching JoelCampo_Burgers.m)
                    if p > 0
                        u_p = u(p);
                    else
                        u_p = conj(u(-p));
                    end
                    
                    if q > 0
                        u_q = u(q);
                    else
                        u_q = conj(u(-q));
                    end
                    
                    conv(j) = conv(j) + u_p * q * 1i * u_q;
                end
            end
            
            % Diffusive term
            if strcmp(case_type, 'DNS')
                % DNS: standard diffusion
                diff(j) = (k(j)^2) * u(j) / params.Re;
            else
                % LES: with eddy viscosity model 
                C_k = case_info.Ck;
                MetLes = 0.31 * (5 - params.m) / (params.m + 1) * sqrt(3 - params.m) * C_k^(-3/2);
                eddy = 1 + 34.5 * exp(-3.03 * (k(N) / k(j)));
                E_k_n = u(N) * conj(u(N));
                nu_eddy = MetLes * (E_k_n / k(N))^0.5 * eddy;
                diff(j) = k(j)^2 * u(j) * (1/params.Re + nu_eddy);
            end
        end
        
        % Update velocities
        u_after = u + time_step * (conv - diff);
        u_after(1) = 1;  % Maintain first mode
        
        % Check convergence
        err = abs(u_after - u);
        if max(err(:)) < params.tolerance
            fprintf("    Converged after %d iterations\n", iter);
            break;
        else
            u = u_after;
        end
        
        if mod(iter, 10000) == 0
            fprintf("    Iteration %d/%d, max error: %.2e\n", iter, params.maxIt, max(err(:)));
        end
    end
    
    if iter == params.maxIt
        fprintf("    Warning: Maximum iterations reached\n");
    end
    
    % Calculate kinetic energy spectrum
    E_k = u_after .* conj(u_after);  % Match JoelCampo_Burgers.m exactly
    
    % Store results
    results(idx).type = case_type;
    results(idx).N = N;
    if strcmp(case_type, 'LES')
        results(idx).Ck = case_info.Ck;
    else
        results(idx).Ck = [];
    end
    results(idx).k = k;
    results(idx).E_k = E_k;
    results(idx).time = iter * time_step;  % Total time to convergence
    results(idx).energy = sum(E_k);
    
    fprintf("  Case %d completed.\n", idx);
end
end
