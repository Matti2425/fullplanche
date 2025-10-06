clearvars

[Geometry,Vin,Pin,Tin,phi_body,gas,NumericalParameters,Init] = readInputData();

[X,Y] = generateMesh(Geometry, NumericalParameters);

S = solidMaterial(Geometry,X,Y);
 

[phi, T, P, rho] = initParams(Init, NumericalParameters, S);


% evaluate coefficients
[a_E, a_W, a_N, a_S, b_P, a_P] = evalCoefficients(phi, phi_body, NumericalParameters,Geometry, Vin, X, Y, S);

% Solver call
[phi, i, error] = solveSystem(phi, phi_body, S, a_E, a_W, a_N, a_S, a_P, b_P, Vin, Y, NumericalParameters, Geometry);

% Compute velocities at faces and nodes
[v_x, v_y, v_mag] = compute_vel(phi, X, Y, NumericalParameters, Vin, Geometry, rho);

% plot(X(S),Y(S),"Color",'red')
% hold on
% plot(X,Y,'o')