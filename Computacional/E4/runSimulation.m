function results = runSimulation(params)
nu = 1 / params.Re;
N_values = params.N_values;
num_cases = numel(N_values);
results = repmat(struct("N", [], "k", [], "E_dns", [], "E_les", [], "time", [], "energy_dns", [], "energy_les", []), num_cases, 1);
for idx = 1:num_cases
    N = N_values(idx);
    fields = initializeSpectralFields(N);
    k_vec = fields.k_vec;
    uk_dns = fields.uk_dns;
    uk_les = fields.uk_les;
    zero_idx = fields.zero_idx;
    mode_ref = fields.mode_reference;
    nu_t = params.Ck * (abs(k_vec) / N) .^ params.m;
    dt = params.CFL / (nu * N^2);
    nt = max(1, round(params.tmax / dt));
    fprintf("  Case %d/%d: N=%d, dt=%.3e, steps=%d\n", idx, num_cases, N, dt, nt);
    time = (0:(nt - 1)).' * dt;
    energy_dns = zeros(nt, 1);
    energy_les = zeros(nt, 1);
    log_interval = max(1, floor(nt / 10));
    for t = 1:nt
        conv_dns = computeConvectiveTerm(k_vec, uk_dns);
        diff_dns = -(k_vec .^ 2) / params.Re .* uk_dns;
        uk_dns = uk_dns + dt * (-conv_dns + diff_dns);
        uk_dns(zero_idx + 1) = mode_ref;
        uk_dns(zero_idx - 1) = conj(mode_ref);
        conv_les = computeConvectiveTerm(k_vec, uk_les);
        diff_les = -(k_vec .^ 2) .* (nu + nu_t) .* uk_les;
        uk_les = uk_les + dt * (-conv_les + diff_les);
        uk_les(zero_idx + 1) = mode_ref;
        uk_les(zero_idx - 1) = conj(mode_ref);
        energy_dns(t) = sum(abs(uk_dns) .^ 2);
        energy_les(t) = sum(abs(uk_les) .^ 2);
        if mod(t, log_interval) == 0 || t == nt
            fprintf("    Progress: %d/%d\n", t, nt);
        end
    end
    k_pos = k_vec((zero_idx + 1):end);
    E_dns = abs(uk_dns((zero_idx + 1):end)) .^ 2;
    E_les = abs(uk_les((zero_idx + 1):end)) .^ 2;
    results(idx).N = N;
    results(idx).k = k_pos;
    results(idx).E_dns = E_dns;
    results(idx).E_les = E_les;
    results(idx).time = time;
    results(idx).energy_dns = energy_dns;
    results(idx).energy_les = energy_les;
    fprintf("  Case N=%d completed.\n", N);
end
end
