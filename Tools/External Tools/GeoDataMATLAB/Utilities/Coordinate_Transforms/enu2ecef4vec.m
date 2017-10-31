function ECEF = enu2ecef4vec(ENU,LatLong)
% enu2ecef4vec.m
% ENU = ecef2enu4vec(ECEF,LatLong)
% by John Swoboda 1/8/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
% This function will take a set of vectors (in terms of the mathematical
% structure) in ENU coordinates and rotate them to an ECEF coordinate
% system.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% ENU - A 3xN matrix (can also be a Nx3 matrix)with the rotated vectors in 
% whatever unit the origial vectors were in.  The east, north and up 
% components are rows 1, 2 and 3 respecitivly.
% LatLon - A 2xN matrix (can also be a Nx2 matrix) with the latitude 
% longitude. The matrix is broken up in the following way, the first row 
% is latitude in degrees, the second row is longitude in degrees.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% ECEF - A 3xN matrix  with X, Y and Z in ECEF coordinate space.  The  
% coordinates can be in what ever units the user wants.  The X,Y and Z 
% coordinates are broken up where they are rows 1, 2, and 3 of the matrix 
% respectively.
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
U = ENU(1,:);
V = ENU(2,:);
W = ENU(3,:);
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
a11 = -sind(long0); a12 = -sind(lat0).*cosd(long0); a13 = cosd(lat0).*cosd(long0);
a21 = cosd(long0); a22 = -sind(lat0).*sind(long0); a23 = cosd(lat0).*sind(long0);
a31 = 0; a32 = cosd(lat0); a33 = sind(lat0);
%% Create vectors

X = a11.*U + a12.*V + a13.*W;
Y = a21.*U + a22.*V + a23.*W;
Z = a31.*U + a32.*V + a33.*W;

ECEF = [X;Y;Z];