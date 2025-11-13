clearvars; close all; clc;
fprintf("Initializing Burgers simulation...\n");
params = setupSimulation();
fprintf("Loaded parameters: Re=%d, cases=%d\n", params.Re, numel(params.cases));

% Try to load cached results
fprintf("Checking for cached results...\n");
[results, cached_params, found] = loadCachedResults(params);

if found
    fprintf("Using cached results (computed previously).\n");
    fprintf("  Cache params: Re=%d, C1=%.6f\n", cached_params.Re, cached_params.C1);
    fprintf("  Current params: Re=%d, C1=%.6f\n", params.Re, params.C1);
    % Verify parameters match (cache filename should ensure this, but double-check)
    if abs(cached_params.Re - params.Re) > 1e-10 || abs(cached_params.C1 - params.C1) > 1e-10
        fprintf("  WARNING: Parameter mismatch! Clearing cache and recomputing...\n");
        cacheFile = getCacheFilename(params);
        if exist(cacheFile, 'file')
            delete(cacheFile);
        end
        found = false;
    else
        % Use cached params in case they differ slightly
        params = cached_params;
    end
end

if ~found
    fprintf("No cache found. Running simulations...\n");
    results = runSimulation(params);
    fprintf("Simulation runs completed. Caching results...\n");
    cacheResults(results, params);
end

fprintf("Generating figures...\n");
plotEnergySpectrum(results, params);
fprintf("All outputs generated.\n");
