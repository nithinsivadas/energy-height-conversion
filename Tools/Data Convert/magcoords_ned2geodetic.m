function [ lat,lon,altitude ] = magcoords_ned2geodetic(magcoords)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Generating lat, lon, h, energy coordinates
xMagEast  = magcoords(:,1);
yMagNorth = magcoords(:,2);
zMagUp    = magcoords(:,3);

lat0 = 65.1260 ;
lon0 = -147.4789;
h0 = 689/1000;  

[lat, lon] = ned2geodetic(yMagNorth, xMagEast, -zMagUp, lat0, lon0, h0, wgs84Ellipsoid('km')); 

altitude = -zMagUp; % altitude and -zMagUp are almost identical. 
end

