function [results, params, found] = loadCachedResults(params, cacheDir)
% LOADCACHEDRESULTS Load simulation results from cache if available
%   [results, params, found] = loadCachedResults(params, cacheDir)
%   attempts to load cached results matching the given parameters.
%
%   Inputs:
%       params  - struct with simulation parameters
%       cacheDir - string path to cache directory (default: 'cache')
%
%   Outputs:
%       results - struct array with simulation results (if found)
%       params  - struct with simulation parameters (from cache)
%       found   - logical, true if cache was found and loaded

if nargin < 2
    cacheDir = 'cache';
end

% Generate cache filename
cacheFile = getCacheFilename(params, cacheDir);

% Check if cache exists
if exist(cacheFile, 'file')
    fprintf('  Loading cached results from: %s\n', cacheFile);
    data = load(cacheFile);
    results = data.results;
    params = data.params;
    found = true;
else
    results = [];
    found = false;
end
end

