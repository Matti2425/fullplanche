clearvars

[Geometry,Vin,Pin,Tin,phi_body,gas,NumericalParameters,Init] = readInputData();

[X,Y] = generateMesh(Geometry, NumericalParameters);

S = solidMaterial(Geometry,X,Y);
 

[phi, T, P, rho] = initParams(Init, NumericalParameters, S, phi_body);

for i = 1:NumericalParameters.N_max
    %evaluate coefficient
    evalCoefficients(phi, phi_body, NumericalParameters,Geometry, Vin, X, Y, S)
    %solve system

    if (X-X_Ans)<error
        break

end



% plot(X(S),Y(S),"Color",'red')
% hold on
% plot(X,Y,'o')