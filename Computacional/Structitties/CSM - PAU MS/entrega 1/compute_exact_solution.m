function [u_exact, du_exact] = compute_exact_solution(A,E,L,c,lambda)

syms x u(x)

b = c*exp(lambda*(x));
Du = diff(u,x);
u_exact = dsolve(diff(u,2)==-b/(E*A), u(0)==0, Du(L)==0);
du_exact = diff(u_exact,x);
end