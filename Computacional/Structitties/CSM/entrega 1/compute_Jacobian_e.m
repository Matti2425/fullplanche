function [Je, dN_dXie] = compute_Jacobian_e(n_G,n_ne,x,e,Tn,Xi_k)

for r=1:1:n_ne
    x_n(r) = x(Tn(e,r));
end
if n_ne == 2
    for jj=1:1:n_G
        for ii=1:1:n_ne
            if ii==1
                dN_dXie(ii,jj) = -1/2;
            else
                dN_dXie(ii,jj) = 1/2;
            end
        end
        Je(1,jj) = dN_dXie(1,jj)*x_n(1) + dN_dXie(2,jj)*x_n(2);
    end
else
    for jj=1:1:n_G
        xi = Xi_k(n_G,jj);
        for ii=1:1:n_ne
            if ii==1
                dN_dXie(ii,jj) = xi-1/2;
            elseif ii==2
                dN_dXie(ii,jj) = -2*xi;
            else
                dN_dXie(ii,jj) = xi+1/2;
            end
        end
        Je(1,jj) = dN_dXie(1,jj)*x_n(1) + dN_dXie(2,jj)*x_n(2) + dN_dXie(3,jj)*x_n(3);
    end
end





end