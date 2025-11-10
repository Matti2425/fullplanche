function sigma_VM = Compute_stress_shells(N_elements,u,Tn,Tm,h_vector,v_poisson_vector,E_vector,Bb,Bmn,Bmt,Bs,R)
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
end

