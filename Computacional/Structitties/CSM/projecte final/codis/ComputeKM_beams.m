function [K,M,f] = ComputeKM_beams(xn,N_nodes,N_elements,Tn,Tm,E_vector,G_vector,rho_vector,j_orientations,A_vector,Iy_vector,Iz_vector,J_vector,ky_vector,kz_vector,kt_vector,Fe,Be,Qe)

    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)

    %--section 2--------------------------------------
    Ndof = 6*N_nodes;
    K=sparse(Ndof,Ndof);
    M=sparse(Ndof,Ndof);

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

end