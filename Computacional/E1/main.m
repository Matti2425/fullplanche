clearvars

[Geometry,Vin,Pin,Tin,phi_body,gas,NumericalParameters,Init] = readInputData();

[X,Y] = generateMesh(Geometry, NumericalParameters);

S = solidMaterial(Geometry,X,Y);
 

[phi, T, P, rho] = initParams(Init, NumericalParameters, S);


% evaluate coefficients
[a_E, a_W, a_N, a_S, b_P, a_P] = evalCoefficients(phi, phi_body, NumericalParameters,Geometry, Vin, X, Y, S, Init);

% Solver call
[phi, i, error] = solveSystem(phi, phi_body, S, a_E, a_W, a_N, a_S, a_P, b_P, Vin, Y, NumericalParameters, Geometry);

% Compute velocities at faces and nodes
[v_x, v_y, v_mag] = compute_vel(phi, X, Y, NumericalParameters, Vin, Geometry, rho);

% Compute thermodynamic properties
[T, P, rho, cp_mean] = compute_thermo(v_mag, NumericalParameters, Init, gas, S);

% Compute aerodynamic coefficients
[C_D, C_L, Gamma, D, L, p_cyl, theta] = final_calc(X, Y, P, v_x, v_y, Geometry, Vin, gas);
figure
%plotCFD('Title',P,NumericalParameters.N,NumericalParameters.M,dx,dy)
plotStreamLines(psi,X,Y)

% plot(X(S),Y(S),"Color",'red')
% hold on
% plot(X,Y,'o')