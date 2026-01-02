function Pe = ComputePressureState(xn,n_u,n_l,b,c,p_inf,alpha)

Pe = sparse(3*(size(n_u,1)+size(n_l,1)),3);

kk = 0;

for ii = 1:size(n_u,1)

    node = n_u(ii,4);
    x = xn(node,1);
    y = xn(node,2);

    if x <= 0.8*b
        Px = p_inf;
    else
        Px = p_inf*cos((pi/2)*(x-0.8*b)/(0.2*b));
    end

    pu = -4*alpha*Px*((1-y/c)^4+sqrt(1-y/c));

    P = -pu*n_u(ii,1:3);
    
    kk = kk+1;
    Pe(kk,:) = [P(1) node 1];
    kk = kk+1;
    Pe(kk,:) = [P(2) node 2];
    kk = kk+1;
    Pe(kk,:) = [P(3) node 3];
end

for ii = 1:size(n_l,1)

    node = n_l(ii,4);
    x = xn(node,1);
    y = xn(node,2);

    if x <= 0.8*b
        Px = p_inf;
    else
        Px = p_inf*cos((pi/2)*(x-0.8*b)/(0.2*b));
    end

    pl = 4*alpha*Px*((1-y/c)^4-0.25*sqrt(1-y/c));

    P = -pl*n_l(ii,1:3);
    
    kk = kk+1;
    Pe(kk,:) = [P(1) node 1];
    kk = kk+1;
    Pe(kk,:) = [P(2) node 2];
    kk = kk+1;
    Pe(kk,:) = [P(3) node 3];
end

end