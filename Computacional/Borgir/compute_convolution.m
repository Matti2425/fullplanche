function C_k = compute_convolution(uk, k_vec, K)

    C_k = zeros(1, K);
    for j = 1:K

        k = k_vec(j);
        sum_triads = 0;

        for p_idx = 1:K

            p = k_vec(p_idx);
            q = k - p;
            q_idx = find(k_vec == q);

            if ~isempty(q_idx)

                up = uk(p_idx);
                uq = uk(q_idx);
                sum_triads = sum_triads + up * (1i * q) * uq;

            end
        end
        C_k(j) = sum_triads;
    end
end