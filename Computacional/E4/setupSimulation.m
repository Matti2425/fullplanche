function params = setupSimulation()
% Define simulation cases:

params.cases = struct('type', {}, 'N', {}, 'Ck', {});
params.cases(1).type = 'DNS';
params.cases(1).N = 20;
params.cases(1).Ck = [];

params.cases(2).type = 'DNS';
params.cases(2).N = 100;
params.cases(2).Ck = [];

params.cases(3).type = 'LES';
params.cases(3).N = 20;
params.cases(3).Ck = 0.4523;

params.cases(4).type = 'LES';
params.cases(4).N = 20;
params.cases(4).Ck = 0.05;

% params.cases(5).type = 'LES';
% params.cases(5).N = 100;
% params.cases(5).Ck = 0.05;
% 
% params.cases(6).type = 'LES';
% params.cases(6).N = 100;
% params.cases(6).Ck = 0.4523;

params.Re = 40;
params.C1 = 0.025;  % Time step coefficient 
params.tolerance = 1e-05;  % Convergence tolerance
params.maxIt = 10^5;  % Maximum iterations
params.m = 2;  % Slope in the energy plot
end
