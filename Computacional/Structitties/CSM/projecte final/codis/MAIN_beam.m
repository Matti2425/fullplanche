%% PART I - BEAM MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all
clear all

%% DATA

problem_type = 1; %1 for bending. 2 for torsion.
Nm = 6; %number of modes

% Define the problem's data (e.g. dimensions, material properties, etc.)
E = 71e9; %Pa
v_poisson = 0.33;
G = E/(2*(1+v_poisson));
rho = 2300; %kg/m3

c = 0.8; %m
b = 5.5; %m
y0 = 0.3083; %m
y1 = 0.176; %m
y2 = 0.464; %m
h1 = 17.5; %m
h2 = 15; %m
h3 = 6; %m
A = 0.0062; %m2
Iy = 8.801e-6; %m4
Iz = 8.069e-5; %m4
J = 8.949e-5; %m4
ky = 0.5457;
kz = 0.1767;
kt = 0.2920;
yc = 0.2954; %m

%% PREPROCESS

% Load mesh data
load('beam.mat','xn','Tn','Tm');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 2]
% Tm : Material connectivities [Nelem x 1]

% Define boundary conditions and forces data matrices: Up, Fe, Qe, Be
N_elements = size(Tn,1);
N_nodes = size(xn,1);
Up = [0 1 1;
      0 1 2;
      0 1 3;
      0 1 4;
      0 1 5;
      0 1 6];

switch(problem_type)
    case 1
    Fe = [1 N_nodes 3];
    case 2
    Fe = [1 N_nodes 4];
end

Qe = [0 1 2];

Be = [0 1 2];

j_orientations = [0 1 0]; 

E_vector = [E]; %more rows would be added if more values existed
A_vector = [A]; %more rows would be added if more values existed
G_vector = [G]; %more rows would be added if more values existed
Iy_vector = [Iy]; %more rows would be added if more values existed
Iz_vector = [Iz]; %more rows would be added if more values existed
ky_vector = [ky]; %more rows would be added if more values existed
kz_vector = [kz]; %more rows would be added if more values existed
kt_vector = [kt]; %more rows would be added if more values existed
J_vector = [J]; %more rows would be added if more values existed
rho_vector = [rho]; %more rows would be added if more values existed


%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
    
    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)

    %--section 2--------------------------------------
    Ndof = 6*N_nodes;
    K=zeros(Ndof,Ndof);
    M=zeros(Ndof,Ndof);

    for e=1:1:N_elements

        l(e) = norm(xn(Tn(e,2),:)-xn(Tn(e,1),:));
        i = ((xn(Tn(e,2),:))' - (xn(Tn(e,1),:))')/l(e);
        j = j_orientations(Tm(e,1),:)';
        k = cross(i,j);
        Rp = [i(1,1) j(1,1) k(1,1) 0 0 0;
             i(2,1) j(2,1) k(2,1) 0 0 0;
             i(3,1) j(3,1) k(3,1) 0 0 0;
             0 0 0 i(1,1) j(1,1) k(1,1);
             0 0 0 i(2,1) j(2,1) k(2,1);
             0 0 0 i(3,1) j(3,1) k(3,1)]';

        R(:,:,e) = [Rp(1,1) Rp(1,2) Rp(1,3) Rp(1,4) Rp(1,5) Rp(1,6) 0 0 0 0 0 0;
                    Rp(2,1) Rp(2,2) Rp(2,3) Rp(2,4) Rp(2,5) Rp(2,6) 0 0 0 0 0 0;
                    Rp(3,1) Rp(3,2) Rp(3,3) Rp(3,4) Rp(3,5) Rp(3,6) 0 0 0 0 0 0;
                    Rp(4,1) Rp(4,2) Rp(4,3) Rp(4,4) Rp(4,5) Rp(4,6) 0 0 0 0 0 0;
                    Rp(5,1) Rp(5,2) Rp(5,3) Rp(5,4) Rp(5,5) Rp(5,6) 0 0 0 0 0 0;
                    Rp(6,1) Rp(6,2) Rp(6,3) Rp(6,4) Rp(6,5) Rp(6,6) 0 0 0 0 0 0;
                    0 0 0 0 0 0 Rp(1,1) Rp(1,2) Rp(1,3) Rp(1,4) Rp(1,5) Rp(1,6);
                    0 0 0 0 0 0 Rp(2,1) Rp(2,2) Rp(2,3) Rp(2,4) Rp(2,5) Rp(2,6);
                    0 0 0 0 0 0 Rp(3,1) Rp(3,2) Rp(3,3) Rp(3,4) Rp(3,5) Rp(3,6);
                    0 0 0 0 0 0 Rp(4,1) Rp(4,2) Rp(4,3) Rp(4,4) Rp(4,5) Rp(4,6);
                    0 0 0 0 0 0 Rp(5,1) Rp(1,2) Rp(5,3) Rp(5,4) Rp(5,5) Rp(5,6);
                    0 0 0 0 0 0 Rp(6,1) Rp(6,2) Rp(6,3) Rp(6,4) Rp(6,5) Rp(6,6)];
        
        Nx(1)=-1/l(e);
        Nx(2)=1/l(e);
        Ba(1,:,e) = [Nx(1) 0 0 0 0 0 Nx(2) 0 0 0 0 0];
        Ca = E_vector(Tm(e,1),1)*A_vector(Tm(e,1),1);
        Ka(:,:,e) = l(e)*R(:,:,e)'*Ba(1,:,e)'*Ca*Ba(1,:,e)*R(:,:,e);

        Bb(:,:,e) = [0 0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0;
                     0 0 0 0 0 Nx(1) 0 0 0 0 0 Nx(2)];

        Cb = E_vector(Tm(e,1),1)*[Iy_vector(Tm(e,1),1) 0;
                                  0 Iz_vector(Tm(e,1),1)];

        Kb(:,:,e) = l(e)*R(:,:,e)'*Bb(:,:,e)'*Cb*Bb(:,:,e)*R(:,:,e);

        N_Gauss = 1/2;

        Bs(:,:,e) = [0 Nx(1) 0 0 0 -N_Gauss 0 Nx(2) 0 0 0 -N_Gauss;
                     0 0 Nx(1) 0 N_Gauss 0 0 0 Nx(2) 0 N_Gauss 0];

        Cs = G_vector(Tm(e,1),1)*A_vector(Tm(e,1),1)*[ky_vector(Tm(e,1),1) 0;
                                                      0 kz_vector(Tm(e,1),1)];

        Ks(:,:,e)=l(e)*R(:,:,e)'*Bs(:,:,e)'*Cs*Bs(:,:,e)*R(:,:,e);

        Bt(1,:,e) = [0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0 0];
        Ct = G_vector(Tm(e,1),1)*J_vector(Tm(e,1),1)*kt_vector(Tm(e,1),1);
        Kt(:,:,e) = l(e)*R(:,:,e)'*Bt(1,:,e)'*Ct*Bt(1,:,e)*R(:,:,e);
        
        Xi = [-1/sqrt(3);
               1/sqrt(3)];
        w = [1;
             1];

        rho_matrix = rho_vector(Tm(e,1),1)*[A_vector(Tm(e,1),1) 0 0 0 0 0;
                                            0 A_vector(Tm(e,1),1) 0 0 0 0;
                                            0 0 A_vector(Tm(e,1),1) 0 0 0;
                                            0 0 0 J_vector(Tm(e,1),1) 0 0;
                                            0 0 0 0 Iy_vector(Tm(e,1),1) 0;
                                            0 0 0 0 0 Iz_vector(Tm(e,1),1)];
        Me(:,:,e)=zeros(12,12);

        for kk=1:1:2
            N(1) = (1-Xi(kk))/2;
            N(2) = (1+Xi(kk))/2;
            N_matrix(:,:,e,kk)=[N(1) 0 0 0 0 0 N(2) 0 0 0 0 0;
                                0 N(1) 0 0 0 0 0 N(2) 0 0 0 0;
                                0 0 N(1) 0 0 0 0 0 N(2) 0 0 0;
                                0 0 0 N(1) 0 0 0 0 0 N(2) 0 0;
                                0 0 0 0 N(1) 0 0 0 0 0 N(2) 0;
                                0 0 0 0 0 N(1) 0 0 0 0 0 N(2)];
            Me(:,:,e) = Me(:,:,e) + (w(kk)*l(e)*R(:,:,e)'*N_matrix(:,:,e,kk)'*rho_matrix*N_matrix(:,:,e,kk)*R(:,:,e))/2;

        end

        for jj=1:1:6
            Idof(jj,1)=6*(Tn(e,1)-1)+jj;
            Idof(6+jj,1)=6*(Tn(e,2)-1)+jj;
        end

        K(Idof,Idof)=K(Idof,Idof) + Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e);
        M(Idof,Idof)=M(Idof,Idof) + Me(:,:,e);

    end


        %--section 3--------------------------------------
        f=zeros(Ndof,1);
        for q=1:1:size(Fe,1)
            f(6*(Fe(q,2)-1)+Fe(q,3),1) = f(6*(Fe(q,2)-1)+Fe(q,3),1) + Fe(q,1);
        end

        Q = zeros(N_nodes,6);

        for r=1:1:size(Qe,1)
            Q(Qe(r,2),Qe(r,3)) = Q(Qe(r,2),Qe(r,3)) + Qe(r,1);
        end

        B = zeros(N_nodes,6);

        for s=1:1:size(Be,1)
            B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3)) + Be(s,1);
        end

        for e=1:1:N_elements
            bb(:,e) = [B(Tn(e,1),:),B(Tn(e,2),:)]';
            qq(:,e) = [Q(Tn(e,1),:),Q(Tn(e,2),:)]';
            fe(:,e) = Me(:,:,e)*bb(:,e);
            for kk=1:1:2
                fe(:,e) = fe(:,e) + (w(kk)*l(e)*R(:,:,e)'*N_matrix(:,:,e,kk)'*N_matrix(:,:,e,kk)*R(:,:,e)*qq(:,e))/2;
            end

            for jj=1:1:6
                Idof(jj,1) = 6*(Tn(e,1)-1)+jj;
                Idof(6+jj,1) = 6*(Tn(e,2)-1)+jj;
            end

            f(Idof,1) = f(Idof,1) + fe(:,e);
        end


        %--section 4--------------------------------------
        u = zeros(Ndof,1);
        for p=1:1:size(Up,1)
            Ip(p) = 6*(Up(p,2)-1) + Up(p,3);
            u(Ip(p),1) = Up(p,1);
        end
        If = setdiff(1:Ndof,Ip);


        %--section 5--------------------------------------
        u(If,1) = K(If,If)\(f(If,1)-(K(If,Ip)*u(Ip,1)));
        fR = K*u-f;

        
        %--section 6--------------------------------------
        for e=1:1:N_elements
            for jj=1:1:6
                Idof(jj,1) = 6*(Tn(e,1)-1)+jj;
                Idof(6+jj,1) = 6*(Tn(e,2)-1)+jj;
            end
            ue = u(Idof,1);

            epsilon_a(1,e) = Ba(1,:,e)*R(:,:,e)*ue;
            epsilon_s(:,e) = Bs(:,:,e)*R(:,:,e)*ue;
            epsilon_t(1,e) = Bt(1,:,e)*R(:,:,e)*ue;
            epsilon_b(:,e) = Bb(:,:,e)*R(:,:,e)*ue;

            fint(:,e) = R(:,:,e)*(Ka(:,:,e)+Kb(:,:,e)+Ks(:,:,e)+Kt(:,:,e))*ue;
            Fx(:,e) = [fint(1,e); -fint(7,e)];
            Fy(:,e) = [fint(2,e); -fint(8,e)];
            Fz(:,e) = [fint(3,e); -fint(9,e)];
            Mx(:,e) = [fint(4,e); -fint(10,e)];
            My(:,e) = [fint(5,e); -fint(11,e)];
            Mz(:,e) = [fint(6,e); -fint(12,e)];
        end


        %--section modal--------------------------------------

    [V,D] = MAIN_modal(K,M,If,Nm);

    Phi = zeros(Ndof,Nm);
    lambda = zeros(1,Nm);
    for kk=1:1:size(V,2)
        Phi(If,kk) = V(:,kk)/sqrt(V(:,kk)'*M(If,If)*V(:,kk));
        for jj=1:1:size(xn,1)
            u_y_modal(jj,kk)=Phi(2+6*(jj-1),kk);
            u_z_modal(jj,kk)=Phi(3+6*(jj-1),kk);
            theta_x_modal(jj,kk)=Phi(4+6*(jj-1),kk);
        end
        lambda(kk) = D(kk,kk);
        freq(kk) = sqrt(lambda(kk))/(2*pi);
    end


    for kk=1:1:size(V,2)
        figure()
        plot(xn(:,1),u_y_modal(:,kk))
        title("u_y mode " + kk)
        figure()
        plot(xn(:,1),u_z_modal(:,kk))
        title("u_z mode " + kk)
        figure()
        plot(xn(:,1),theta_x_modal(:,kk))
        title("theta_x mode " + kk)
    end

    
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

for i=1:1:N_nodes
    u_z(i,1) = u(3+6*(i-1),1);
end
for i=1:1:N_nodes
    theta_x(i,1) = u(4+6*(i-1),1);
end


figure()
plot(xn,u_z)
title('Vertical deflection along the spanwise')
switch(problem_type)
    case 1
        legend('bending test at beam model','Location','northwest')
    case 2
        legend('torsion test at beam model','Location','northwest')
end
xlabel('x (m)')
ylabel('u_z (m)')

figure()
plot(xn,theta_x)
title('Twist angle along the spanwise')
switch(problem_type)
    case 1
        legend('bending test at beam model','Location','northwest')
    case 2
        legend('torsion test at beam model','Location','northwest')
end
xlabel('x (m)')
ylabel('\theta_x (rad)')











