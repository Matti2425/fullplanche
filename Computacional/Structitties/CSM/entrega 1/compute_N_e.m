function N = compute_N_e(Xi_k,n_G,i,k,n_ne)

if n_ne == 2
    if i==1
        N = (1/2)*(1-Xi_k(n_G,k));
    elseif i==2
        N = (1/2)*(1+Xi_k(n_G,k));
    end
else
    if i==1
        N = (1/2)*Xi_k(n_G,k)*(Xi_k(n_G,k)-1);
    elseif i==2
        N = 1-(Xi_k(n_G,k))^2;
    elseif i==3
        N = (1/2)*Xi_k(n_G,k)*(Xi_k(n_G,k)+1);
    end
end



end