function plotEnergyEvolution(results)
fprintf("  Plotting total energy evolution (%d cases)...\n", numel(results));
figure
hold on
for idx = 1:numel(results)
    data = results(idx);
    plot(data.time, data.energy_dns, "-", "LineWidth", 1.5, "DisplayName", sprintf("DNS N=%d", data.N))
    plot(data.time, data.energy_les, "--", "LineWidth", 1.5, "DisplayName", sprintf("LES N=%d", data.N))
end
xlabel("Time")
ylabel("Total energy")
title("Total energy evolution")
grid on
legend("Location", "best")
hold off
fprintf("  Energy evolution plot complete.\n");
end
