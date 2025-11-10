function conv = computeConvectiveTerm(k_vec, uk)
K = numel(k_vec);
conv = zeros(K, 1);
for j = 1:K
    k = k_vec(j);
    sum_triads = 0;
    for p_idx = 1:K
        p = k_vec(p_idx);
        q = k - p;
        q_idx = find(k_vec == q, 1);
        if ~isempty(q_idx)
            sum_triads = sum_triads + uk(p_idx) * 1i * q * uk(q_idx);
        end
    end
    conv(j) = sum_triads;
end
end
