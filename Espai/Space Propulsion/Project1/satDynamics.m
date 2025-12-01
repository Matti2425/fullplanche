function dy = satDynamics(t,y,mu,F,ce,alpha_max)
    r=y(1); th=y(2); ph=y(3);
    Vr=y(4); Vt=y(5); Vp=y(6);
    m=y(7);

    V = sqrt(Vr^2+Vt^2+Vp^2);
    alpha = alpha_max*cos(ph);

    % F_in = F*cos(alpha);
    % F_perp = F*sin(alpha);
    % 
    % Fr = F_in*(Vr/max(V,1e-12));
    % % Ftheta = 0;
    % Fphi = F_perp;

    Ftheta = F*cos(alpha);
    Fr = F*sin(alpha);
    Fphi = 0; %?  
   
    dy=zeros(7,1);

    dy(1)=Vr;
    dy(2)=Vt/r;
    dy(3)=Vp/(r*sin(th));

    dy(4)= -mu/r^2 + (Vt^2+Vp^2)/r + Fr/m;
    dy(5)= -(Vr*Vt)/r + Vp^2/(r*tan(th)) + Ftheta/m;
    dy(6)= -(Vr*Vp)/r - (Vt*Vp)/(r*tan(th)) + Fphi/m;

    dy(7)= -F/ce;
end