function dy = satDynamics(t,y,mu,F,ce,alpha_max)
    r=y(1); th=y(2); ph=y(3);
    Vr=y(4); Vt=y(5); Vp=y(6);
    m=y(7);

    alpha = alpha_max*cos(ph);

    % F_in = F*cos(alpha);
    % F_perp = F*sin(alpha);
    % 
    % Fr = F_in*(Vr/max(V,1e-12));
    % % Ftheta = 0;
    % Fphi = F_perp;

    Ftheta = F*cos(alpha);
    Fr = 0;
    Fphi = F*sin(alpha); %?  
   

    dy = [Vr;
        Vt/r;
        Vp/(r*sin(th));
        -mu/r^2 + (Vt^2+Vp^2)/r + Fr/m;
        -(Vr*Vt)/r + Vp^2/(r*tan(th)) + Ftheta/m;
        -(Vr*Vp)/r - (Vt*Vp)/(r*tan(th)) + Fphi/m;
        -F/ce];
end