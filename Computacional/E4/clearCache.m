function clearCache(cacheDir)
% CLEARCACHE Clear all cached simulation results
%   clearCache(cacheDir) deletes all .mat files in the cache directory.
%   If cacheDir is not specified, defaults to 'cache'.
%
%   Example:
%       clearCache()  % Clear default cache directory
%       clearCache('my_cache')  % Clear custom cache directory

if nargin < 1
    cacheDir = 'cache';
end

if ~exist(cacheDir, 'dir')
    fprintf('Cache directory "%s" does not exist. Nothing to clear.\n', cacheDir);
    return;
end

% Find all .mat files in cache directory
cacheFiles = dir(fullfile(cacheDir, '*.mat'));

if isempty(cacheFiles)
    fprintf('No cached files found in "%s".\n', cacheDir);
    return;
end

% Delete all cache files
deleted = 0;
for i = 1:length(cacheFiles)
    filePath = fullfile(cacheDir, cacheFiles(i).name);
    delete(filePath);
    deleted = deleted + 1;
end

fprintf('Cleared %d cached file(s) from "%s".\n', deleted, cacheDir);
end

