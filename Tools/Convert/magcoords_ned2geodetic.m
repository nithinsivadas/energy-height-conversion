function [ lat,lon,altitude ] = magcoords_ned2geodetic(magcoords, lat0,lon0, h0)
%% magcoords_ned2geodetic Convert magnetic field aligned ned coordinates to geodetic at PFISR
%--------------------------------------------------------------------------
% Input
%------
% magcoords - [nCoordinates x nDimenson]
%     (:,1) - East  [km]
%     (:,2) - North [km]
%     (:,3) - Up    [km]
%--------------------------------------------------------------------------
% Output
%-------
% lat      - latitude in deg [nCoordinates x 1]
% lon      - longtude in deg [nCoordinates x 1]
% altitude - in km [nCoordinates x 1]
%--------------------------------------------------------------------------
% Modified: 2nd Feb 2018 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------

%% Generating lat, lon, h, energy coordinates
xMagEast  = magcoords(:,1);
yMagNorth = magcoords(:,2);
zMagUp    = magcoords(:,3);

if nargin<4
    h0 = 689/1000;  
end
if nargin<3
    lon0 = -147.4789;
end
if nargin<2
    lat0 = 65.1260 ;
end


[lat, lon] = ned2geodetic(yMagNorth, xMagEast, -zMagUp, lat0, lon0, h0, wgs84Ellipsoid('km')); 

altitude = -zMagUp; % altitude and -zMagUp are almost identical. 
end

