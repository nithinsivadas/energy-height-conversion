function ECEF = enu2ecefl(ENU,LatLongHeight)
% enu2ecef.m
% ENU = ecef2enu(ECEF,LatLongHeight)
% by John Swoboda 1/8/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
% This function will take a set of east north up vectors with the origin
% points in LatLongHeight.  LatLongHeight can be a single set of
% coordinates or a collection of coordintes to deal with a moving platform.
% This function uses enu2ecef4vec.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% ENU - A 3xN matrix (can also be a Nx3 matrix)with the rotated vectors in 
% meters from the origin point stated in LatLongHeight.  The east, north 
% and up components are rows 1, 2 and 3 respecitivly.
% LatLonHeight - A 3xN matrix (can also be a Nx3 matrix) with the latitude 
% longitude and height. These will be the coordinates that the ENU 
% coordinate system is in reference to.  The matrix is broken up in the 
% following way, the first row is latitude in degrees, the second row is 
% longitude in degrees and the last row is the height in meters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% ECEF - A 3xN matrix with X, Y and Z in ECEF coordinate space.  The  
% coordinates must be in meters.  The X,Y and Z coordinates are broken up 
% where they are rows 1, 2, and 3 of the matrix respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
% Wikipedia article: http://en.wikipedia.org/wiki/Geodetic_system
%% Data check
% Reorient data so its in a 1xN arrays 
[~,n] = size(ENU);
data_trans = false;

if (n ==3)
    ENU = ENU';
    data_trans=true;
end

% assume that the lat and long data are oriented the same way
if data_trans == true
    LatLongHeight = LatLongHeight';
end
%% Perfom Calculation
ECEF0 = wgs2ecef(LatLongHeight);
LatLong = LatLongHeight(1:2,:);
% correct for a single point
if size(LatLongHeight,2)==1
    ECEF0=repmat(ECEF0,1,size(ENU,2));
end
ECEF = enu2ecef4vec(ENU,LatLong)+ECEF0;
