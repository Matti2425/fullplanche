clearvars

[Geometry,vin,pin,Tin,phi_body,gas,NumericalParameters] = readInputData();

[X,Y] = generateMesh(Geometry, NumericalParameters);

evalCoefficients(System)

S = solidMaterial(Geometry,X,Y);



% plot(X(S),Y(S),"Color",'red')
% hold on
% plot(X,Y,'o')