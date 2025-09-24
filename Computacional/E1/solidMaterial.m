function S = solidMaterial(Geometry,X,Y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    dist = sqrt((X - Geometry.Center(1)).^2 + (Y - Geometry.Center(2)).^2);

    % S = 1 si dintre S = 0 si fora
    S = dist < Geometry.D/2;

end