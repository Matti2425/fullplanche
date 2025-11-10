%% PART II - WING MODELLING 

% Useful commands to initialize script
clc
clear
close all

%% SOLUTION CONTROLS

problem_type = 3; % 1 for bending, 2 for torsion, 3 for pressure and weight
solution_type = 2; % 1 for direct method, 2 for model-order reduction solution
Nm = 9; % Number of desired modes to be computed (this is the number that will be used to solve the system if model-order reduction solution is requested)
Im = [1,2,3,4,6,7]; %Number of modes to be used for the modal projection

WB = 1; %Set to 1 to include, set to 0 to neglect
ST = 1; %Set to 1 to include, set to 0 to neglect
RB = 1; %Set to 1 to include, set to 0 to neglect
SK = 1; %Set to 1 to include, set to 0 to neglect

computeModes = WB*ST*RB*SK; % Don't compute modes if any is 0

if computeModes == 0 && solution_type == 2
    error('Model-order reduction is only implemented for the solving of the entire model, it does not work if parts are neglected.')
end


%% DATA

% Define the problem's data (e.g. dimensions, material properties, etc.)
E = 71e9; %Pa
v_poisson = 0.33; %chan chachachan chachachchanchanchan chan chachachachan chachac
rho = 2300; %kg/m3
G = E/(2*(1+v_poisson));

c = 0.8; %m
b = 5.5; %m
y0 = 0.3083; %m
y1 = 0.176; %m
y2 = 0.464; %m
yc = 0.2954; %m

% Loads 
F = 1; %N
T = 1; %Nm

% other

p_inf = 12500; % Pa
alpha = 10*pi/180; % rad
g0 = 9.80665; %m/s^2



%% PREPROCESS

% Load mesh data
load('wing.mat','xn','Tn_st','Tm_st','Tn_wb','Tm_wb','Tn_rb','Tm_rb','Tn_sk','Tm_sk','indRoot','indPointA','indPointB','indSpar1','indSpar2','n_u','n_l');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn_st : Nodal connectivities for beam elements (st: stringers) [Nelem x 2]
% Tn_wb, Tn_rb, Tn_sk : Nodal connectivities for shell elements (wb: wingbox, rb: ribs, sk: skin) [Nelem x 4]
% Tm_st, Tm_wb, Tm_rb, Tm_sk : Material connectivities for the different elements [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.
% n_u, n_l  : Matrices containing information about the unit normals in the upper 
%             and lower surfaces, respectively. 
%       |nx1 ny1 nz1 id1|   First three columns are components of normal vector 
%       |nx2 ny2 nz2 id2|   Last (fourth) column is the nodal index
%       |      ...      |

% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be

N_nodes = size(xn,1);
Ndof = 6*N_nodes;

% LOAD PREPROCESSING

kk=0;
for i=1:1:size(indRoot,1)
    for j=1:1:6
        kk=kk+1;
        Up(kk,1) = 0;
        Up(kk,2) = indRoot(i,1);
        Up(kk,3) = j;
    end
end

Fa_shear = F*((y2-yc)/(y2-y1));
Fb_shear = F*((yc-y1)/(y2-y1));
Fa_torque = -T/(y2-y1);
Fb_torque = T/(y2-y1);


% WING BOX

N_elements_wb = size(Tn_wb,1);

h1_wb = 17.5e-3; %m
h2_wb = 15e-3; %m
h3_wb = 6e-3; %m
h_vector_wb = [h1_wb; h2_wb; h3_wb]; %more rows would be added if more values existed

E_vector_wb = [E;E;E];%more rows would be added if more values existed
v_poisson_vector_wb = [v_poisson; v_poisson; v_poisson]; %more rows would be added if more values existed
rho_vector_wb = [rho; rho; rho]; %more rows would be added if more values existed

switch(problem_type)
    case 1
        Fe_wb = [Fa_shear indPointA 3;
              Fb_shear indPointB 3];
        Pe_wb = [];
    case 2
        Fe_wb = [Fa_torque indPointA 3;
              Fb_torque indPointB 3];
        Pe_wb = [];
    case 3
        Fe_wb = [];
        Pe_wb = [];

end



% STRINGERS

N_elements_st = size(Tn_st,1);

ky_st = 5/6;
kz_st = 5/6;
kt_st = 1;

r_st = 5*10^-3; %m

A_st = pi*(r_st)^2;%more rows would be added if more values existed


Iy_st = (pi/4)*r_st^4; %m^4
Iz_st = (pi/4)*r_st^4; %m^4

J_st = (pi/2)*r_st^4; %m^4

E_vector_st = E;%more rows would be added if more values existed

G_vector_st = G;%more rows would be added if more values existed

rho_vector_st = rho;%more rows would be added if more values existed

A_vector_st = A_st; %more rows would be added if more values existed

Iy_vector_st = Iy_st;
Iz_vector_st = Iz_st;
J_vector_st = J_st;

ky_vector_st = ky_st;
kz_vector_st = kz_st;
kt_vector_st = kt_st;

j_orientations = [0 1 0]; 

switch(problem_type)
    case 1
        Fe_st = [];
        Qe_st = [];
    case 2
        Fe_st = [];
        Qe_st = [];
    case 3
        Fe_st = [];
        Qe_st = [];
end

% RIBS

N_elements_rb = size(Tn_rb,1);

h1_rb = 5e-3; %m

h_vector_rb = [h1_rb]; %more rows would be added if more values existed

E_vector_rb = [E];%more rows would be added if more values existed
v_poisson_vector_rb = [v_poisson]; %more rows would be added if more values existed
rho_vector_rb = [rho]; %more rows would be added if more values existed

switch(problem_type)
    case 1
        Fe_rb = [];
        Pe_rb = [];
    case 2
        Fe_rb = [];
        Pe_rb = [];
    case 3
        Fe_rb = [];
        Pe_rb = [];
end

% SKIN

N_elements_sk = size(Tn_sk,1);

h1_sk = 2e-3; %m

h_vector_sk = [h1_sk]; %more rows would be added if more values existed

E_vector_sk = [E];%more rows would be added if more values existed
v_poisson_vector_sk = [v_poisson]; %more rows would be added if more values existed
rho_vector_sk = [rho]; %more rows would be added if more values existed

switch(problem_type)
    case 1
        Fe_sk = [];
        Pe_sk = [];
    case 2
        Fe_sk = [];
        Pe_sk = [];
    case 3
        Fe_sk = [];
        Pe_sk = ComputePressureState(xn,n_u,n_l,b,c,p_inf,alpha);
end

% BODY FORCES
switch(problem_type)
    case 1
        Be = [];
    case 2
        Be = [];
    case 3
        for ii = 1:size(xn,1)

              Be(ii,:) = [-g0 ii 3];

        end

end


%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if  1% Set to 1 to (re)compute the system matrices and '0' to load them
    
    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)
    % ...
    [K_wb,M_wb,f_wb,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb,R_wb] = ComputeKM_shells(xn,N_nodes,N_elements_wb,Tn_wb,Tm_wb,h_vector_wb,E_vector_wb,v_poisson_vector_wb,rho_vector_wb,Fe_wb,Be,Pe_wb);
    [K_st,M_st,f_st] = ComputeKM_beams(xn,N_nodes,N_elements_st,Tn_st,Tm_st,E_vector_st,G_vector_st,rho_vector_st,j_orientations,A_vector_st,Iy_vector_st,Iz_vector_st,J_vector_st,ky_vector_st,kz_vector_st,kt_vector_st,Fe_st,Be,Qe_st);
    [K_rb,M_rb,f_rb,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb,R_rb] = ComputeKM_shells(xn,N_nodes,N_elements_rb,Tn_rb,Tm_rb,h_vector_rb,E_vector_rb,v_poisson_vector_rb,rho_vector_rb,Fe_rb,Be,Pe_rb);
    [K_sk,M_sk,f_sk,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk,R_sk] = ComputeKM_shells(xn,N_nodes,N_elements_sk,Tn_sk,Tm_sk,h_vector_sk,E_vector_sk,v_poisson_vector_sk,rho_vector_sk,Fe_sk,Be,Pe_sk);

    % Once (re)computed, save them to a separate data file
    save('wing_matrices.mat','K_wb','M_wb','K_st','M_st','K_rb','M_rb','K_sk','M_sk','f_wb','f_st','f_rb','f_sk');
    save('wing_aux_matrices.mat','Bb_wb','Bmn_wb','Bmt_wb','Bs_wb','R_wb','Bb_rb','Bmn_rb','Bmt_rb','Bs_rb','R_rb','Bb_sk','Bmn_sk','Bmt_sk','Bs_sk','R_sk');
    % TIP: Add other potential results that can be reused in other parts
    % (e.g. element's length 'l', elements rotations matrices 'R', etc.)
else
    
    % Load previously computed results
    load('wing_matrices.mat','K_wb','M_wb','K_st','M_st','K_rb','M_rb','K_sk','M_sk','f_wb','f_st','f_rb','f_sk');
    load('wing_aux_matrices.mat','Bb_wb','Bmn_wb','Bmt_wb','Bs_wb','R_wb','Bb_rb','Bmn_rb','Bmt_rb','Bs_rb','R_rb','Bb_sk','Bmn_sk','Bmt_sk','Bs_sk','R_sk');
end


% Compute global matrices

K = WB*K_wb + ST*K_st + RB*K_rb + SK*K_sk;
M = WB*M_wb + ST*M_st + RB*M_rb + SK*M_sk;

% Compute global force vector

f = f_wb + f_st + f_rb + f_sk;

% Boundary conditions

    u = zeros(Ndof,1);
    for pp=1:1:size(Up,1)
        Ip(pp)=6*(Up(pp,2)-1)+Up(pp,3);
        u(Ip(pp),1) = Up(pp,1);
    end
    If=setdiff(1:Ndof,Ip);

% Perform modal analysis

if computeModes == 1
[V,D] = MAIN_modal(K,M,If,Nm);

    Phi = zeros(Ndof,Nm);
    lambda = zeros(1,Nm);
    for kk=1:1:size(V,2)
        Phi(If,kk) = V(:,kk)/sqrt(V(:,kk)'*M(If,If)*V(:,kk));
        lambda(kk) = D(kk,kk);
        freq(kk) = sqrt(lambda(kk))/(2*pi);
    end
end

% Solve system

if solution_type == 1 % Direct solution
    u(If,1)=K(If,If)\(f(If,1)-K(If,Ip)*u(Ip,1));
    fR = K*u-f;
elseif solution_type == 2 % Model-order reduction solution
    F = f; % more columns would be added if different frequencies were analyzed. In this case, only wk = 0 is analyzed (static).
    omega = 0; % more columns would be added if different frequencies were analyzed. In this case, only wk = 0 is analyzed (static).
    u = zeros(Ndof,length(omega));
    for kk = 1:size(F,2)
        for jj = 1:length(Im)
            alpha_modal(jj,kk) = (Phi(:,Im(jj)))'*F(:,kk)/(lambda(jj)-omega(kk)^2);
            u(:,kk) = u(:,kk) + Phi(:,Im(jj))*alpha_modal(jj,kk);
        end
    end

else
    error('Requested solution type does not exist. The currently implemented solution types are: 1 for direct solving; 2 for solving with model-order reduction.')
end

% Posterior calculations

sigVM_wb = Compute_stress_shells(N_elements_wb,u,Tn_wb,Tm_wb,h_vector_wb,v_poisson_vector_wb,E_vector_wb,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb,R_wb);
sigVM_rb = Compute_stress_shells(N_elements_rb,u,Tn_rb,Tm_rb,h_vector_rb,v_poisson_vector_rb,E_vector_rb,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb,R_rb);
sigVM_sk = Compute_stress_shells(N_elements_sk,u,Tn_sk,Tm_sk,h_vector_sk,v_poisson_vector_sk,E_vector_sk,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk,R_sk);




%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
%save('wing_results.mat');

% Include plot functions

for i=1:1:size(indSpar1,1)
    theta_x(i,1) = (u(3+6*(indSpar2(i,1)-1),1) - u(3+6*(indSpar1(i,1)-1),1))/(y2-y1);
    u_z(i,1) = u(3+6*(indSpar1(i,1)-1),1) + theta_x(i,1)*(yc-y1);
    xn_plot(i,1) = xn(indSpar1(i,1),1);
end


figure()
plot(xn_plot,u_z)
title('Vertical deflection along the spanwise')
switch(problem_type)
    case 1
        legend('bending test at wing model','Location','northwest')
    case 2
        legend('torsion test at wing model','Location','northwest')
    case 3
        legend('wing model under pressure and body forces','Location','northwest')
end
xlabel('x (m)')
ylabel('u_z (m)')

figure()
plot(xn_plot,theta_x)
title('Twist angle along the spanwise')
switch(problem_type)
    case 1
        legend('bending test at wing model','Location','northwest')
    case 2
        legend('torsion test at wing model','Location','northwest')
    case 3
        legend('wing model under pressure and body forces','Location','northwest')
end
xlabel('x (m)')
ylabel('\theta_x (rad)')

% Additional plot functions useful to visualize 3D model and modes
scale = 1;
plotDeformed('wing',xn,Tn_wb,u,scale,sigVM_wb); % For wingbox elements
plotDeformed('wing',xn,Tn_rb,u,scale,sigVM_rb); % For rib elements
plotDeformed('wing',xn,Tn_sk,u,scale,sigVM_sk); % For skin elements

% This function plots the deformed structure and Von Mises stress distribution: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).
% sigVM : Von Mises stress at each Gauss point. It is expected to be given as 
%         a matrix with dimensions [Nelem x Ngauss].
imodes = [1,2,3,4,5,6,7,8,9];

plotModes('wing',Phi,freq,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 