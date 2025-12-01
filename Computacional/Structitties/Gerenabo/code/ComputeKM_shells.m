function [K,M,f,Bb,Bmn,Bmt,Bs,R] = ComputeKM_shells(xn,N_nodes,N_elements,Tn,Tm,h_vector,E_vector,v_poisson_vector,rho_vector,Fe,Be,Pe)
   
%--section 2--------------------------------------
    Ndof = 6*N_nodes;
    K=sparse(Ndof,Ndof);
    M=sparse(Ndof,Ndof);

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


% Compute external forces vector

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
end