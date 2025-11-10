function fields = initializeSpectralFields(N)
k_vec = (-N:N).';
K = numel(k_vec);
zero_idx = N + 1;
uk = zeros(K, 1);
for j = 1:N
    value = 1 / j;
    uk(zero_idx + j) = value;
    uk(zero_idx - j) = conj(value);
end
uk(zero_idx) = 0;
fields.k_vec = k_vec;
fields.uk_dns = uk;
fields.uk_les = uk;
fields.zero_idx = zero_idx;
fields.mode_reference = uk(zero_idx + 1);
end
