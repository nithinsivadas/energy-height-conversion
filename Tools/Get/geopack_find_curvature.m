function [Kc] = geopack_find_curvature(XX,YY,ZZ)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

dX = gradient(XX);
ddX = gradient(dX);

dY = gradient(YY);
ddY = gradient(dY);

dZ = gradient(ZZ);
ddZ = gradient(dZ);

Kc = sqrt((ddZ.*dY - ddY.*dZ).^2 + (ddX.*dZ - ddZ.*dX).^2 + (ddY.*dX - ddX.*dY).^2)...
.*(dX.^2 + dY.^2 + dZ.^2).^-1.5;

end
