function filename = getCacheFilename(params, cacheDir)
% GETCACHEFILENAME Generate a unique cache filename based on parameters
%   filename = getCacheFilename(params, cacheDir) creates a filename
%   based on the case parameters to uniquely identify the simulation.
%
%   Inputs:
%       params  - struct with simulation parameters
%       cacheDir - string path to cache directory
%
%   Output:
%       filename - full path to cache file

if nargin < 2
    cacheDir = 'cache';
end

% Create a string representation of all cases
caseStr = '';
for idx = 1:numel(params.cases)
    case_info = params.cases(idx);
    if strcmp(case_info.type, 'LES')
        caseStr = sprintf('%s_%s_N%d_Ck%.6f', caseStr, case_info.type, ...
            case_info.N, case_info.Ck);
    else
        caseStr = sprintf('%s_%s_N%d', caseStr, case_info.type, case_info.N);
    end
end

% Include other important parameters
paramStr = sprintf('Re%d_C1%.6f_tol%.2e_maxIt%d_m%d', ...
    params.Re, params.C1, params.tolerance, params.maxIt, params.m);

% Create hash-like filename (remove leading underscore if present)
if ~isempty(caseStr) && caseStr(1) == '_'
    caseStr = caseStr(2:end);
end
fullStr = [caseStr '_' paramStr];

% Replace problematic characters for filename
fullStr = regexprep(fullStr, '[\.\+\-\s]', '_');
fullStr = regexprep(fullStr, '_+', '_');  % Replace multiple underscores with single

filename = fullfile(cacheDir, [fullStr '.mat']);
end

