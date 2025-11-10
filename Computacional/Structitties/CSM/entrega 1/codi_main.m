%
close all
clear
clc


%-------input data---------
% % A = 120; %mm2
% % E = 75; %GPa
% % L = 2; %m
% % c = 1e4; %N/m
A = 120*10^(-6); %m2
E = 75*10^9; %Pa
L = 2; %m
c = 1e4; %N/m
lambda = 2;

% n_ne = 2;
%n_G = 2;


w_k = [2 0 0;
       1 1 0;
       5/9 8/9 5/9];

Xi_k = [0 0 0;
       -1/sqrt(3) 1/sqrt(3) 0;
       -sqrt(3/5) 0 sqrt(3/5)];

%exact solutions
[u_exact, du_exact] = compute_exact_solution(A,E,L,c,lambda);
x_vector_exact = linspace(0,L,500);
u_vector_exact = (x_vector_exact*exp(4))/1800 - exp(2*x_vector_exact)/3600 + 1/3600;
sigma_vector_exact = E*(exp(4)/1800 - exp(2*x_vector_exact)/1800);
sigma_vector_exact=sigma_vector_exact*10^(-6); %conversion to Mpa


%finite element solution
for sol=1:1:2
    
    clear u

    if sol == 1
        n_G = 1;
        n_ne = 2;
    elseif sol == 2
        n_G = 2;
        n_ne = 3;
    end

for cont_n_el = 1:1:7
    n_el = 2^(cont_n_el);
    n = n_ne + (n_el-1)*(n_ne-1);

    h = L/n_el;
    h_vector(cont_n_el) = h;

    x = linspace(0,L,n)';

    Tn = zeros(n_el,n_ne);

    for e=1:1:n_el
        for i = 1:1:n_ne
            Tn(e,i) = (e-1)*(n_ne-1) + i;
        end
    end

    B = zeros(1,n_ne,n_el,n_G);
    Ke = zeros(n_ne,n_ne,n_el);
    K = zeros(n,n);

    
    for e=1:1:n_el
        le = x(Tn(e,n_ne))-x(Tn(e,1));
        [Je, dN_dXie] = compute_Jacobian_e(n_G,n_ne,x,e,Tn,Xi_k);
        for k=1:1:n_G
            for i=1:1:n_ne
                B(1,i,e,k)=(1/Je(k))*dN_dXie(i,k);
            end
            Ke(:,:,e)=Ke(:,:,e)+w_k(n_G,k)*Je(k)*E*A*(B(:,:,e,k)'*B(:,:,e,k));
        end
        for i=1:1:n_ne
            for j=1:1:n_ne
                K(Tn(e,i),Tn(e,j)) = K(Tn(e,i),Tn(e,j)) + Ke(i,j,e);
            end
        end
    end


    N = zeros(1,n_ne,n_el,n_G);
    x_G = zeros(n_G,n_el);
    fe = zeros(n_ne,n_el);
    f = zeros(n,1);

    for e=1:1:n_el
        xe = x(Tn(e,:));
        le = x(Tn(e,n_ne))-x(Tn(e,1));
        [Je, dN_dXie] = compute_Jacobian_e(n_G,n_ne,x,e,Tn,Xi_k);
        for k=1:1:n_G
            for i=1:1:n_ne
                N_i = compute_N_e(Xi_k,n_G,i,k,n_ne);
                N(1,i,e,k)=N_i;
            end

            x_G(k,e) = N(:,:,e,k)*xe;
            b_k = c*exp(lambda*(x_G(k,e)));
            fe(:,e) = fe(:,e) + w_k(n_G,k)*Je(k)*N(:,:,e,k)'*b_k;
        end
        for i=1:1:n_ne
            f(Tn(e,i),1) = f(Tn(e,i),1) + fe(i,e);
        end
    end


    v_r = 1;
    v_f = setdiff(1:n,v_r);

    K_ff = K(v_f,v_f);
    K_fr = K(v_f,v_r);
    K_rf = K(v_r,v_f);
    K_rr = K(v_r,v_r);

    u(v_r,1) = 0;
    u(v_f,1) = K_ff\(f(v_f,1)-K_fr*u(v_r,1));
    f_r = K_rr*u(v_r,1) + K_rf*u(v_f,1)-f(v_r,1);



    sigma = zeros(n_G,n_el);
    sigma_G=zeros(n_G,n_el);
    for e=1:1:n_el
        u_e = u(Tn(e,:),1);
        for k=1:1:n_G
            sigma_G(k,e) = E*B(:,:,e,k)*u_e;
        end
    end

    sigma=interp1(x_G(:),sigma_G(:),x,'linear','extrap');

    sigma=sigma*10^(-6); %Mpa
if sol == 1
    figure()
    plot(x,u,':o');
    hold on;
    plot(x_vector_exact,u_vector_exact)
    title("Displacements along the 1D beam for " + n_el + " linear elements")
    xlabel('x (m)')
    ylabel('u (m)')
    legend('finite element solution','exact solution','Location','northeast')
    hold off;
    
    figure()
    plot(x,sigma,':o');
    hold on;
    plot(x_vector_exact,sigma_vector_exact)
    title("Stress along the 1D beam for " + n_el + " linear elements")
    xlabel('x (m)')
    ylabel('\sigma (MPa)')
    legend('finite element solution','exact solution','Location','northeast')
    hold off;

elseif sol == 2
    figure()
    plot(x,u,':o');
    hold on;
    plot(x_vector_exact,u_vector_exact)
    title("Displacements along the 1D beam for " + n_el + " quadratic elements")
    xlabel('x (m)')
    ylabel('u (m)')
    legend('finite element solution','exact solution','Location','northeast')
    hold off;
    figure()
    plot(x,sigma,':o');
    hold on;
    plot(x_vector_exact,sigma_vector_exact)
    title("Stress along the 1D beam for " + n_el + " quadratic elements")
    xlabel('x (m)')
    ylabel('\sigma (MPa)')
    legend('finite element solution','exact solution','Location','northeast')
    hold off;

end

u_exact_error=(x.*exp(4))/1800 - exp(2.*x)/3600 + 1/3600;

e_numerador = sqrt(sum((u_exact_error-u).^2));
e_denominador = sqrt(sum((u_exact_error).^2));
error(cont_n_el,n_G) = e_numerador/e_denominador;


end


end


figure()
loglog(h_vector,error(:,1),'-o')
hold on
loglog(h_vector,error(:,2),'-.+')
title("Accuracy plot")
xlabel('h (m)')
ylabel('error')
legend('Linear elements','Quadratic elements','Location','southeast')







