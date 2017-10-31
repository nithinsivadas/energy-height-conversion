function ENU = ecef2enul(ECEF,LatLongHeight)
% ecef2enul.m
% ENU = ecef2enul(ECEF,LatLongHeight)
% by John Swoboda 1/10/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
% This fuction will take a set of ecef coorinates in meters and transfers
% them to an east north up coordinate system also in meters.  LatLongHeight 
% can be a single set of coordinates or a collection of coordintes to deal 
% with a moving platform.This function uses ecef2enu4vec.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% ECEF - A 3xN matrix (can also be a Nx3 matrix) with X, Y and Z in ECEF 
% coordinate space.  The  coordinates must be in meters.  The X,Y and Z 
% coordinates are broken up where they are rows 1, 2,
% and 3 of the matrix respectively(assuming a 3xN array, if transposed this
% will be the columns.
% LatLonHeight - A 3xN matrix (can also be a Nx3 matrix) with the latitude 
% longitude and height. These will be the coordinates that the ENU 
% coordinate system is in reference to.  The matrix is broken up in the 
% following way, the first row is latitude in degrees, the second row is 
% longitude in degrees and the last row is the height in meters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% ENU - A 3xN matrix with the rotated vectors in meters refereced to east 
% up from the original stated in the WGS coordinate given.  The east, north
% and up components are rows 1, 2 and 3 respectivly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
% Wikipedia article: http://en.wikipedia.org/wiki/Geodetic_system
%% Data check
% Reorient data so its in a 1xN arrays 
[~,n] = size(ECEF);
data_trans = false;

if (n ==3)
    ECEF = ECEF';
    data_trans=true;
end

% assume that the lat and long data are oriented the same way
if data_trans == true
    LatLongHeight = LatLongHeight';
end
%% Perfom Calculation
ECEF0 = wgs2ecef(LatLongHeight);
LatLong = LatLongHeight(1:2,:);
ENU = ecef2enu4vec(ECEF-ECEF0,LatLong);
