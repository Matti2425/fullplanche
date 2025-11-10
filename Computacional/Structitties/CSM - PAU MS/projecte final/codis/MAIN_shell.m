%% PART I - SHELL MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all
clear all

%% DATA

problem_type = 2; %1 for bending. 2 for torsion.
Nm = 6; %number of modes

% Define the problem's data (e.g. dimensions, material properties, etc.)
E = 71e9; %Pa
v_poisson = 0.33;
G = E/(2*(1+v_poisson));
rho = 2300; %kg/m3

F = 1; %N
T = 1; %Nm

c = 0.8; %m
b = 5.5; %m
y0 = 0.3083; %m
y1 = 0.176; %m
y2 = 0.464; %m
h1 = 17.5e-3; %m
h2 = 15e-3; %m
h3 = 6e-3; %m
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
load('shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 4]
% Tm : Material connectivities [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.

% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be
N_elements = size(Tn,1);
N_nodes = size(xn,1);

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

switch(problem_type)
    case 1
        Fe = [Fa_shear indPointA 3;
              Fb_shear indPointB 3];
    case 2
        Fe = [Fa_torque indPointA 3;
              Fb_torque indPointB 3];
end

Pe = [0 1 2];
Be = [0 1 2];

h_vector = [h1; h2; h3]; %more rows would be added if more values existed
E_vector = [E; E; E]; %more rows would be added if more values existed
v_poisson_vector = [v_poisson; v_poisson; v_poisson]; %more rows would be added if more values existed
rho_vector = [rho; rho; rho]; %more rows would be added if more values existed



%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
    
    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)


    %--section 2--------------------------------------
    Ndof = 6*N_nodes;
    K=sparse(Ndof,Ndof);
    M=sparse(Ndof,Ndof);
    %K=zeros(Ndof,Ndof);
    %M=zeros(Ndof,Ndof);

    for e=1:1:N_elements
        term1=xn(Tn(e,3),:)'-xn(Tn(e,1),:)';
        term2=xn(Tn(e,4),:)'-xn(Tn(e,2),:)';
        S = cross(term1,term2)/2;
        k_dir = S/norm(S);
        d = (xn(Tn(e,2),:)'+xn(Tn(e,3),:)'-xn(Tn(e,4),:)'-xn(Tn(e,1),:)')/2;
        i_dir = d/norm(d);
        j_dir = cross(k_dir,i_dir);
        Rp=[i_dir(1,1) j_dir(1,1) k_dir(1,1) 0 0;
            i_dir(2,1) j_dir(2,1) k_dir(2,1) 0 0;
            i_dir(3,1) j_dir(3,1) k_dir(3,1) 0 0;
            0 0 0 i_dir(1,1) j_dir(1,1);
            0 0 0 i_dir(2,1) j_dir(2,1);
            0 0 0 i_dir(3,1) j_dir(3,1)]';

        R(:,:,e)=zeros(20,24);
        R(1:5,1:6,e)=Rp;
        R(6:10,7:12,e)=Rp;
        R(11:15,13:18,e)=Rp;
        R(16:20,19:24,e)=Rp;

        a_coef = [-1 1 1 -1];
        b_coef = [-1 -1 1 1];
        N1=[1 1 1 1]'/4;
        N1_xi = a_coef/4;
        N1_nu = b_coef/4;
        J1 = zeros(2,2);

        for i=1:1:4
            matrix_Ni = [N1_xi(i); N1_nu(i)];
            matrix_ij = [i_dir j_dir];
            J1 = J1 + matrix_Ni*xn(Tn(e,i),:)*matrix_ij;
        end
        N1x = J1\[N1_xi; N1_nu];
        S1=4*det(J1);

        for i=1:1:4
            Bs_i(:,:,i)=[0 0 N1x(1,i) 0 N1(i);
                       0 0 N1x(2,i) -N1(i) 0];
        end

        Cs = eye(2)*5*h_vector(Tm(e,1),1)*E_vector(Tm(e,1),1)/(12*(1+v_poisson_vector(Tm(e,1),1)));
        Bs(:,:,e) = [Bs_i(:,:,1) Bs_i(:,:,2) Bs_i(:,:,3) Bs_i(:,:,4)];
        Ks(:,:,e) = S1*R(:,:,e)'*Bs(:,:,e)'*Cs*Bs(:,:,e)*R(:,:,e);

        for i=1:1:4
            Bmt_i(:,:,i) = [N1x(2,i) N1x(1,i) 0 0 0];
        end
        Cmt = h_vector(Tm(e,1),1)*E_vector(Tm(e,1),1)/(2*(1+v_poisson_vector(Tm(e,1),1)));
        Bmt(:,:,e) = [Bmt_i(:,:,1) Bmt_i(:,:,2) Bmt_i(:,:,3) Bmt_i(:,:,4)];
        Km(:,:,e) = S1*R(:,:,e)'*Bmt(:,:,e)'*Cmt*Bmt(:,:,e)*R(:,:,e);

        Kb(:,:,e) = zeros(24,24);
        Me(:,:,e) = zeros(24,24);
        xi4 = [-1 1 1 -1]/sqrt(3);
        nu4 = [-1 -1 1 1]/sqrt(3);
        w4 = [1 1 1 1];

        for kk=1:1:4
            J4 = zeros(2,2);
            for i=1:1:4
                N4(i)=(1+a_coef(i)*xi4(kk))*(1+b_coef(i)*nu4(kk))/4;
                N4_xi(1,i)=a_coef(i)*(1+b_coef(i)*nu4(kk))/4;
                N4_nu(1,i)=b_coef(i)*(1+a_coef(i)*xi4(kk))/4;
                matrix_N4 = [N4_xi(i); N4_nu(i)];
                J4 = J4 + matrix_N4*xn(Tn(e,i),:)*matrix_ij;
            end
            N4x = J4\[N4_xi; N4_nu];
            S4(e,kk) = w4(kk)*det(J4);
            for i=1:1:4
                Bmn_i(:,:,i)=[N4x(1,i) 0 0 0 0;
                            0 N4x(2,i) 0 0 0];
            end
            Cmn = [1 v_poisson_vector(Tm(e,1),1); v_poisson_vector(Tm(e,1),1) 1]*h_vector(Tm(e,1),1)*E_vector(Tm(e,1),1)/(1-v_poisson_vector(Tm(e,1),1)^2);
            Bmn(:,:,e,kk) = [Bmn_i(:,:,1) Bmn_i(:,:,2) Bmn_i(:,:,3) Bmn_i(:,:,4)];
            Km(:,:,e) = Km(:,:,e) + S4(e,kk)*R(:,:,e)'*Bmn(:,:,e,kk)'*Cmn*Bmn(:,:,e,kk)*R(:,:,e);
            for i=1:1:4
                Bb_i(:,:,i) = [0 0 0 0 N4x(1,i);
                             0 0 0 N4x(2,i) 0;
                             0 0 0 -N4x(1,i) N4x(2,i)];
            end
            Cb = [1 v_poisson_vector(Tm(e,1),1) 0; v_poisson_vector(Tm(e,1),1) 1 0; 0 0 (1-v_poisson_vector(Tm(e,1),1))/2]*h_vector(Tm(e,1),1)^3*E_vector(Tm(e,1),1)/(12*(1-v_poisson_vector(Tm(e,1),1)^2));
            Bb(:,:,e,kk) = [Bb_i(:,:,1) Bb_i(:,:,2) Bb_i(:,:,3) Bb_i(:,:,4)];
            Kb(:,:,e) = Kb(:,:,e) + S4(e,kk)*R(:,:,e)'*Bb(:,:,e,kk)'*Cb*Bb(:,:,e,kk)*R(:,:,e);
            for i=1:1:4
                N_i(:,:,i) = N4(i)*eye(5);
            end
            rho_matrix = rho_vector(Tm(e,1),1)*h_vector(Tm(e,1),1)*[1 0 0 0 0;
                                                                    0 1 0 0 0;
                                                                    0 0 1 0 0;
                                                                    0 0 0 h_vector(Tm(e,1),1)^2/12 0;
                                                                    0 0 0 0 h_vector(Tm(e,1),1)^2/12];
            N(:,:,e,kk) = [N_i(:,:,1) N_i(:,:,2) N_i(:,:,3) N_i(:,:,4)];
            Me(:,:,e) = Me(:,:,e) + S4(e,kk)*R(:,:,e)'*N(:,:,e,kk)'*rho_matrix*N(:,:,e,kk)*R(:,:,e);
        end
        for jj=1:1:6
            Idof(jj,1) = 6*(Tn(e,1)-1)+jj;
            Idof(6+jj,1) = 6*(Tn(e,2)-1)+jj;
            Idof(12+jj,1) = 6*(Tn(e,3)-1)+jj;
            Idof(18+jj,1) = 6*(Tn(e,4)-1)+jj;
        end
        K(Idof,Idof) = K(Idof,Idof) + Km(:,:,e) + Kb(:,:,e) + Ks(:,:,e);
        M(Idof,Idof) = M(Idof,Idof) + Me(:,:,e);

    end

    %--section 3--------------------------------------
    n = zeros(3,N_nodes);
    for e=1:1:N_elements
        term11=xn(Tn(e,3),:)'-xn(Tn(e,1),:)';
        term22=xn(Tn(e,4),:)'-xn(Tn(e,2),:)';
        SS = cross(term11,term22)/2;
        SS_mod(e) = sqrt(SS(1)^2+SS(2)^2+SS(3)^2);
        kp(:,e) = SS/SS_mod(e);
        for i=1:1:size(Tn,2)
            n(:,Tn(e,i)) = n(:,Tn(e,i)) + kp(:,e);
        end
    end

    %Kr=zeros(Ndof,Ndof);
    Kr=sparse(Ndof,Ndof);
    for e=1:1:N_elements
        for i=1:1:size(Tn,2)
            alpha = acos(dot(n(:,Tn(e,i)),kp(:,e))/norm(n(:,Tn(e,i))));
            alpha_deg = alpha*360/(2*pi);
            if alpha_deg<0.1
                Idof = 6*(Tn(e,i)-1)+[4 5 6]';
                Kr(Idof,Idof) = Kr(Idof,Idof) + E_vector(Tm(e,1),1)*h_vector(Tm(e,1),1)*SS_mod(e)*kp(:,e)*kp(:,e)';
            end
        end
    end

    K=K+Kr;

    
    %--section 4--------------------------------------
    f = zeros(Ndof,1);
    for qq=1:1:size(Fe,1)
        f(6*(Fe(qq,2)-1)+Fe(qq,3),1)=f(6*(Fe(qq,2)-1)+Fe(qq,3),1) + Fe(qq,1);
    end
    P=zeros(N_nodes,6);
    for rr=1:1:size(Pe,1)
        P(Pe(rr,2), Pe(rr,3))=P(Pe(rr,2), Pe(rr,3)) + Pe(rr,1);
    end
    B=zeros(N_nodes,6);
    for sss=1:1:size(Be,1)
        B(Be(sss,2), Be(sss,3))=B(Be(sss,2), Be(sss,3)) + Be(sss,1);
    end

    for e=1:1:N_elements
        bbb(:,e) = [B(Tn(e,1),:) B(Tn(e,2),:) B(Tn(e,3),:) B(Tn(e,4),:)]';
        p(:,e) = [P(Tn(e,1),:) P(Tn(e,2),:) P(Tn(e,3),:) P(Tn(e,4),:)]';
        fe(:,e) = Me(:,:,e)*bbb(:,e);
        for kk=1:1:4
            fe(:,e) = fe(:,e) + S4(e,kk)*R(:,:,e)'*N(:,:,e,kk)'*N(:,:,e,kk)*R(:,:,e)*p(:,e);
        end
        for jj=1:1:6
            Idof(jj,1) = 6*(Tn(e,1)-1)+jj;
            Idof(6+jj,1) = 6*(Tn(e,2)-1)+jj;
            Idof(12+jj,1) = 6*(Tn(e,3)-1)+jj;
            Idof(18+jj,1) = 6*(Tn(e,4)-1)+jj;
        end
        f(Idof,1)=f(Idof,1)+fe(:,e);
    end


    %--section 5--------------------------------------
    u = zeros(Ndof,1);
    for pp=1:1:size(Up,1)
        Ip(pp)=6*(Up(pp,2)-1)+Up(pp,3);
        u(Ip(pp),1) = Up(pp,1);
    end
    If=setdiff(1:Ndof,Ip);


    %--section 6--------------------------------------
    u(If,1)=K(If,If)\(f(If,1)-K(If,Ip)*u(Ip,1));
    fR = K*u-f;


    %--section 7--------------------------------------
    for e=1:1:N_elements
        for jj=1:1:6
            Idof(jj,1) = 6*(Tn(e,1)-1)+jj;
            Idof(6+jj,1) = 6*(Tn(e,2)-1)+jj;
            Idof(12+jj,1) = 6*(Tn(e,3)-1)+jj;
            Idof(18+jj,1) = 6*(Tn(e,4)-1)+jj;
        end
        for kk=1:1:4
            epsilon_b(:,e,kk)=Bb(:,:,e,kk)*R(:,:,e)*u(Idof,1);
            epsilon_m(1:2,e,kk)=Bmn(:,:,e,kk)*R(:,:,e)*u(Idof,1);
            epsilon_m(3,e,kk)=Bmt(:,:,e)*R(:,:,e)*u(Idof,1);
            epsilon_s(:,e,kk)=Bs(:,:,e)*R(:,:,e)*u(Idof,1);
        end
        Cp=[1 v_poisson_vector(Tm(e,1),1) 0; v_poisson_vector(Tm(e,1),1) 1 0; 0 0 (1-v_poisson_vector(Tm(e,1),1))/2]*E_vector(Tm(e,1),1)/(1-v_poisson_vector(Tm(e,1),1)^2);
        Cs=eye(2)*E_vector(Tm(e,1),1)/(2*(1+v_poisson_vector(Tm(e,1),1)));
        for kk=1:1:4
            sigma_m(:,e,kk)=Cp*epsilon_m(:,e,kk);
            sigma_s(:,e,kk)=Cs*epsilon_s(:,e,kk);
            sigma_b(:,e,kk)=Cp*h_vector(Tm(e,1),1)*epsilon_b(:,e,kk)/2;
            sigma_plus=[sigma_m(:,e,kk)+sigma_b(:,e,kk); sigma_s(:,e,kk)]';
            sigma_VM_plus=(sigma_plus(1)^2+sigma_plus(2)^2-sigma_plus(1)*sigma_plus(2)+3*(sigma_plus(3)+sigma_plus(4)+sigma_plus(5)))^(1/2);
            sigma_minus=[sigma_m(:,e,kk)-sigma_b(:,e,kk); sigma_s(:,e,kk)]';
            sigma_VM_minus=(sigma_minus(1)^2+sigma_minus(2)^2-sigma_minus(1)*sigma_minus(2)+3*(sigma_minus(3)+sigma_minus(4)+sigma_minus(5)))^(1/2);
            sigma_VM(e,kk)=max(sigma_VM_plus,sigma_VM_minus);
        end
    end

            %--section modal--------------------------------------

    [V,D] = MAIN_modal(K,M,If,Nm);

    Phi = zeros(Ndof,Nm);
    lambda = zeros(1,Nm);
    for kk=1:1:size(V,2)
        Phi(If,kk) = V(:,kk)/sqrt(V(:,kk)'*M(If,If)*V(:,kk));
        for jj=1:1:size(indSpar1,1)
            u_y_modal(jj,kk) = (Phi(2+6*(indSpar1(jj,1)-1),kk)+Phi(2+6*(indSpar2(jj,1)-1),kk))/2;
            theta_x_modal(jj,kk) = (Phi(3+6*(indSpar2(jj,1)-1),kk) - Phi(3+6*(indSpar1(jj,1)-1),kk))/(y2-y1);
            u_z_modal(jj,kk) = Phi(3+6*(indSpar1(jj,1)-1),kk) + theta_x_modal(jj,kk)*(yc-y1);
        end
        lambda(kk) = D(kk,kk);
        freq(kk) = sqrt(lambda(kk))/(2*pi);
    end

        for kk=1:1:size(V,2)
        figure()
        plot(xn(indSpar1(:),1),u_y_modal(:,kk))
        title("u_y mode " + kk)
        figure()
        plot(xn(indSpar1(:),1),u_z_modal(:,kk))
        title("u_z mode " + kk)
        figure()
        plot(xn(indSpar1(:),1),theta_x_modal(:,kk))
        title("\theta_x mode " + kk)
        end    


    
    % Once (re)computed, save them to a separate data file
    save('shell_matrices.mat','K','M'); 
    % TIP: Add other potential results that can be reused in other parts
    % (e.g. element's length 'l', elements rotations matrices 'R', etc.)
    
else
    
    % Load previously computed results
    load('shell_matrices.mat','K','M');
    
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
save('shell_results.mat');

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
        legend('bending test at shell model','Location','northwest')
    case 2
        legend('torsion test at shell model','Location','northwest')
end
xlabel('x (m)')
ylabel('u_z (m)')

figure()
plot(xn_plot,theta_x)
title('Twist angle along the spanwise')
switch(problem_type)
    case 1
        legend('bending test at shell model','Location','northwest')
    case 2
        legend('torsion test at shell model','Location','northwest')
end
xlabel('x (m)')
ylabel('\theta_x (rad)')














% Additional plot functions useful to visualize 3D model and modes
%scale = 1;
% plotDeformed('shell',xn,Tn,u,scale);
% This function plots the deformed structure: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).
%imodes = [1,2,3,4,5,6];
% plotModes('shell',Phi,freq,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 