% ANTONIO LUQUE LÓPEZ
% MUEA 
% ASSIGNMENT 1: POTENTIAL FLOW IN A CHANNEL
%--------------------------------------------------------------------------
clear workspace,
clear all;
clc;        
clearvars;  
close all;

%--------------------------------------------------------------------------

%%INSTRUCTIONS%%

%CASE=1 for channel; CASE=2 for channel+static cylinder; CASE=3 for channel+rotating cylinder(w>0-->Clockwise)
%SLINES=1 Plot streamlines(ONLY FOR v_P, v_xP, v_yP); SLINES=0 Don't plot streamlines
%MESH=1 -->Plot the mesh(NOT RECOMENDED); MESH=0 -->Don't plot the mesh
%PLOT=1 for v_P, 2 for v_xP, 3 for v_yP, 4 for p_P, 5 for T_P-->Plot the variable field in a colormap.

CASE=2; SLINES=1; MESH=0; PLOT=1;

%--------------------------------------------------------------------------

% 1. INPUT DATA
% 1.1 Physical data
L=2;
H=1;
v_in=1;
p_in=101325;
T_in=300;

gamma = 1.4;   
R = 287;

% 1.2 Numerical data

N = 5;       
M = 7;        
delta = 1e-10;
psi_start = 0; 

% 1.2 Reference data

p_ref=101325;
T_ref=300;
v_ref=1;
rho_F = p_in / (R * T_in);
rho_ref=rho_F;

% 1.3 Cylinder data

if CASE==3
    w=25; % Enter an angular velocity
else
    w=0;
end

if CASE~=1
    R_cyl=0.1;
    x_center=L/2;
    y_center=H/2;
    theta=0;
    psi_body = v_in*(H/2)+R_cyl^2*w*(1/(2*pi));% Valid for static and rotating cylinder
else
    R_cyl=0;
    x_cyl=0;
    x_center=0;
    y_center=0;
    psi_body =0;
end

% 2. PREVIOUS CALCULATIONS

dx = L/N;   
dy = H/M;   
[x_P, y_P] = mesh_generator(N,M,L,H); % Function that creates the geometry of the mesh.

% 3.INITIAL ESTIMATED VALUES
psi = psi_start * ones(M+2, N+2);
solid = ones(M+2, N+2); %1 for fluid, 0 for solid.

if CASE~=1
 [x_cyl, y_cyl, solid, theta] = cylinder_generator(x_center, y_center, R_cyl, M, N, x_P, y_P, solid); %Function that creates the cylinder.
end

rho = rho_F*solid; %Density rho_F for fluid, 0 for solid.

% 4. EVALUATE THE DISCRETIZATION COEFFICIENTS
[a_E, a_W, a_N, a_S, b_P, a_P] = compute_coef(psi_start, psi_body, M, N, L, H, v_in, y_P, x_P, solid);

% 5. SOLVE THE SET OF DISCRETIZED EQUATIONS
[psi, iter, error] = gauss_seidel(psi, psi_body, solid, a_E, a_W, a_N, a_S, a_P, b_P, v_in, y_P, M, N, H, delta);

% 6. OBTAIN THE VELOCITIES AT THE FACES AND THE NODES
[v_xP, v_yP, v_P] = compute_velocity(psi, y_P, x_P, M, N, v_in, rho_ref, rho);

% 7. COMPUTE THERMODYNAMICS PROPERTIES
[T_P, p_P, rho_P, cp_mean] = compute_thermo_properties(T_ref, p_ref, v_ref, v_P, R, gamma, M, N);

% 8. FINAL CALCULATIONS
if CASE~=1
    [C_D, C_L, Gamma, D, L, p_cyl] = final_calculations(x_P, y_P, p_P, v_xP, v_yP, x_cyl, y_cyl, theta, R_cyl, w, rho_ref, v_in);
end
% 9. PRINT RESULTS

switch PLOT
    case 1
        field = v_P;
        title = 'v_P';
    case 2
        fieldToPlot = v_xP;
        title = 'v_xP';
    case 3
        field = v_yP;
        title = 'v_yP';
    case 4
        field = p_P;
        title = 'p_P';
    case 5
        field = T_P;
        title = 'T_P';
end

plotCFD(title, field, N, M, dx, dy);
hold on;

if SLINES==1 % Print the streamlines.
plotStreamLines(psi, x_P, y_P) 
hold on
end

if MESH==1 % Print the mesh.
mesh_plot(x_P, y_P, L, H, N, M) 
hold on
end

if CASE~=1 %Print the cylinder
    fill(x_center + (R_cyl+0.7*dx) * cos(linspace(0, 2*pi, 100)),y_center + (R_cyl+0.7*dx) * sin(linspace(0, 2*pi, 100)), 'k', 'FaceAlpha', 1, 'EdgeColor', 'none'); 
end


