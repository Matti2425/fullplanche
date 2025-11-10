function [V,D] = MAIN_modal(K,M,If,Nm)
K = (K+K')/2;
M = (M+M')/2;
[V,D] = eigs(K(If,If),M(If,If),Nm,'sm');
end