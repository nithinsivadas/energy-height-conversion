function ENU = ecef2enu4vec(ECEF,LatLong)
% ecef2enu4vec.m
% ENU = ecef2enu4vec(ECEF,LatLong)
% by John Swoboda 1/8/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
% This function will take a set of vectors (in terms of the mathematical
% structure) in ECEF coordinates and rotate them to an ENU coordinate
% system.  This compares well with the internal matlab function and yeilds
% a normalized error close to 10^-12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% ECEF - A 3xN matrix (can also be a Nx3 matrix) with X, Y and Z in ECEF 
% coordinate space.  The  coordinates can be in what ever units the user 
% wants.  The X,Y and Z coordinates are broken up where they are rows 1, 2,
% and 3 of the matrix respectively(assuming a 3xN array, if transposed this
% will be the columns.
% LatLon - A 2xN matrix (can also be a Nx2 matrix) with the latitude 
% longitude. The matrix is broken up in the following way, the first row 
% is latitude in degrees, the second row is longitude in degrees.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% ENU - A 3xN matrix with the rotated vectors in whatever unit the origial
% vectors were in.  The east, north and up components are rows 1, 2 and 3
% respecitivly.
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
U = ECEF(1,:);
V = ECEF(2,:);
W = ECEF(3,:);
% assume that the lat and long data are oriented the same way
if data_trans == true
    LatLong = LatLong';
end
% Get the latitude and longitude into 
if numel(LatLong)==2
    lat0 = LatLong(1);
    long0 = LatLong(2);
    lat0 = lat0(ones(1,length(U)));
    long0 = long0(ones(1,length(U)));
else
    lat0 = LatLong(1,:);
    long0 = LatLong(2,:);
end


%% Set up calculation
a11 = -sind(long0); a12 = cosd(long0); a13 = zeros(size(lat0));
a21 = -sind(lat0).*cosd(long0); a22 = -sind(lat0).*sind(long0); a23 = cosd(lat0);
a31 = cosd(lat0).*cosd(long0); a32 = cosd(lat0).*sind(long0); a33 = sind(lat0);
%% Create vectors

East = a11.*U + a12.*V + a13.*W;
North = a21.*U + a22.*V + a23.*W;
Up = a31.*U + a32.*V + a33.*W;

ENU = [East;North;Up];