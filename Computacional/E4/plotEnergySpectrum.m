function plotEnergySpectrum(results, params)
fprintf("  Plotting energy spectrum (%d cases)...\n", numel(results));
figure
hold on
set(gca, "XScale", "log", "YScale", "log")
set(gca, "LineWidth", 1)

% Plot each case
for idx = 1:numel(results)
    case_info = results(idx);
    if strcmp(case_info.type, 'DNS')
        loglog(case_info.k, case_info.E_k, "o-", "LineWidth", 1.5, ...
            "DisplayName", sprintf("DNS N=%d", case_info.N))
    else
        loglog(case_info.k, case_info.E_k, "s-", "LineWidth", 1.5, ...
            "DisplayName", sprintf("LES N=%d Ck=%.4f", case_info.N, case_info.Ck))
    end
end

% Add reference slope - use the case with maximum N to span full range
max_N = 0;
max_N_idx = 1;
for idx = 1:numel(results)
    if results(idx).N > max_N
        max_N = results(idx).N;
        max_N_idx = idx;
    end
end

ref_k = results(1).k(1);
ref_val = results(1).E_k(1);
k_ref = results(max_N_idx).k;  % Use k values from case with maximum N
ref_curve = (k_ref / ref_k) .^ (-params.m) * ref_val;
loglog(k_ref, ref_curve, "k--", "LineWidth", 1.2, ...
    "DisplayName", sprintf("Slope -%d", params.m))

xlabel("Mode k", 'FontSize', 16)
ylabel("E_k", 'FontSize', 16)
title(sprintf("Energy spectrum, Re=%d", params.Re), "FontSize", 16)
legend("Location", "southwest", "FontSize", 16)
grid on
hold off
fprintf("  Energy spectrum plot complete.\n");
end
