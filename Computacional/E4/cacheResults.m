function cacheResults(results, params, cacheDir)
% CACHERESULTS Save simulation results to cache directory
%   cacheResults(results, params, cacheDir) saves the results structure
%   and params to a cache file based on a hash of the case parameters.
%
%   Inputs:
%       results - struct array with simulation results
%       params  - struct with simulation parameters
%       cacheDir - string path to cache directory (default: 'cache')

if nargin < 3
    cacheDir = 'cache';
end

% Create cache directory if it doesn't exist
if ~exist(cacheDir, 'dir')
    mkdir(cacheDir);
end

% Generate cache filename based on parameters
cacheFile = getCacheFilename(params, cacheDir);

% Save results and params
save(cacheFile, 'results', 'params', '-v7.3');
fprintf('  Results cached to: %s\n', cacheFile);
end

