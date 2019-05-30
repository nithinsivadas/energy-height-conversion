function [distance] = geodetic_distance(coord1,coord2,radius)
%% geodetic_distance.m Shortest distance between two points in geodetic coordinates
%           using Harversine formulae 
%   Input:
%   coord1 & coord2 - lat and lon (East longitude)
%   radius - [m] distance from the center to the point
% https://www.movable-type.co.uk/scripts/latlong.html
%
%   Output:
%   distance in [m]

%% Calculation using Harversine Formulae

coordRad1 = deg2rad(coord1);
coordRad2 = deg2rad(coord2);
lambda1 = coordRad1(:,1); phi1 = coordRad1(:,2);
lambda2 = coordRad2(:,1); phi2 = coordRad2(:,2);
dphi = phi2-phi1;
dlambda = lambda2-lambda1;
a = (sin(dphi./2)).^2 + cos(phi1).*cos(phi2).*(sin(dlambda./2)).^2;
c = 2.*atan2(a.^0.5,(1-a).^0.5);
distance = radius.*c;


end

