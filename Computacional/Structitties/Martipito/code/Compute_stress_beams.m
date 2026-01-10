function sigVM_st = Compute_stress_beams(xn, Tn_st, u, E, G, r, kt)
% sigVM_st: [Nelem x 1] von Mises-like measure at outer fiber for each stringer element

Nelem = size(Tn_st,1);
sigVM_st = zeros(Nelem,1);

for e = 1:Nelem
    n1 = Tn_st(e,1);
    n2 = Tn_st(e,2);

    x1 = xn(n1,:)'; 
    x2 = xn(n2,:)';
    Le = norm(x2 - x1);

    Nx = [-1/Le, 1/Le];

    % Extract nodal DOFs for the element (only what we need)
    ux  = [u(6*(n1-1)+1); u(6*(n2-1)+1)];
    thx = [u(6*(n1-1)+4); u(6*(n2-1)+4)];
    thy = [u(6*(n1-1)+5); u(6*(n2-1)+5)];
    thz = [u(6*(n1-1)+6); u(6*(n2-1)+6)];

    eps_x   = Nx * ux;          % dux/dx
    dthx_dx = Nx * thx;         % d(theta_x)/dx
    ky      = Nx * thy;         % d(theta_y)/dx  (curvature proxy)
    kz      = Nx * thz;         % d(theta_z)/dx  (curvature proxy)

    % Outer-fiber normal stress (axial + bending)
    sigma_max = E * (eps_x + r*sqrt(ky^2 + kz^2));

    % Outer-fiber torsional shear (using torsion stiffness definition T = G*J*kt*dthx/dx)
    tau_max = G * kt * r * dthx_dx;

    % Von Mises
    sigVM_st(e) = sqrt(sigma_max^2 + 3*tau_max^2);
end
end