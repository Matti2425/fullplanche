%% PART I - BEAM MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all

%% DATA

% Define the problem's data (e.g. dimensions, material properties, etc.)
% ...

%% PREPROCESS
zeros(0,3);
% Load mesh data
load('beam.mat','xn','Tn','Tm');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 2]
% Tm : Material connectivities [Nelem x 1]

% Define boundary conditions and forces data matrices: Up, Fe, Qe, Be
% ...

%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
    
    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)
    % ...
    
    % Once (re)computed, save them to a separate data file
    save('beam_matrices.mat','K','M'); 
    % TIP: Add other potential results that can be reused in other parts
    % (e.g. element's length 'l', elements rotations matrices 'R', etc.)
    
else
    
    % Load previously computed results
    load('beam_matrices.mat','K','M');
    
end

% Perform modal analysis
% ...

% Compute external forces vector
% ...

% Solve system
% ...

%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('beam_results.mat');

% Include plot functions
% ...