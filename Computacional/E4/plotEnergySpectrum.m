function plotEnergySpectrum(results, params)
fprintf("  Plotting energy spectrum (%d cases)...\n", numel(results));
figure
hold on
set(gca, "XScale", "log", "YScale", "log")
set(gca, "LineWidth", 1)
for idx = 1:numel(results)
    loglog(results(idx).k, results(idx).E_dns, "o-", "LineWidth", 1.5, "DisplayName", sprintf("DNS N=%d", results(idx).N))
    loglog(results(idx).k, results(idx).E_les, "s-", "LineWidth", 1.5, "DisplayName", sprintf("LES N=%d", results(idx).N))
end
ref_k = results(1).k(1);
ref_val = results(1).E_dns(1);
ref_curve = (results(end).k / ref_k) .^ (-params.m) * ref_val;
loglog(results(end).k, ref_curve, "k--", "LineWidth", 1.2, "DisplayName", sprintf("Slope -%d", params.m))
xlabel("Mode k")
ylabel("E_k")
title(sprintf("DNS vs LES spectra, Re=%d", params.Re))
legend("Location", "southwest")
grid on
hold off
fprintf("  Energy spectrum plot complete.\n");
end
