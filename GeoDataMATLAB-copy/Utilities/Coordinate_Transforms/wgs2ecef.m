function ECEF_COORDS = wgs2ecef(WGS_COORDS)
% wgs2ecef.m
% ECEF_COORDS = wgs2ecef(WGS_COORDS)
% by John Swoboda 1/8/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
% This function takes a set of latitude longitude height coordinates and
% puts them into ECEF.  The input is assumed to be a 3xN (or Nx3) matrix of
% WGS coordinates.  It then outputs a 3xN matrix of ECEF coordinates in
% meters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% WGS_COORDS - A 3xN matrix (can also be a Nx3 matrix) with the latitude
% longitude and height coordinates.  The matrix is broken up in the
% following way, the first row is latitude in degrees, the second row is
% longitude in degrees and the last row is height from the surface of the
% earth in meters.   This has been compared to matlab's internal
% geodetic2ecef and has comparible results with a normalized difference on 
% the order of 10^-16. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output 
% ECEF_COORDS -A 3xN matrix with X, Y and Z in ECEF coordinate space.  The 
% coordinates are in units of meters. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference
% Vermeille H (2002) Direct transformation from geocentric coordinates
% to geodetic coordinates, Journal of Geodesy, vol. 76, no. 8, pp 451 - 454,
% Nov 2002, http://link.springer.com/article/10.1007%2Fs00190-002-0273-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input
[n,m] = size(WGS_COORDS);
data_trans = false;

if n~=3&&m==3
    WGS_COORDS = WGS_COORDS.';
    data_trans = true;

elseif n==3&&m==3
    warning('Both dimensions have 3 elements, assuming lat, long and height are along rows');
elseif n~=3 &&m~=3
    error('Neither dimension has a size of 3');
end
%% Get Lat, Long and Height
phi = WGS_COORDS(1,:);
lambda = WGS_COORDS(2,:);
h = WGS_COORDS(3,:);
%% Set the constants
a = 6378137; % semi-major axis in meters
f = 1/298.257223563; % the flattening factor
b = a*(1-f);% semiminor axis in meters 

e = sqrt(a^2-b^2)/a;% first eccentricity

M_e = (a*(1-e^2))./(1-e^2*sind(phi)).^(3/2);% Meridian radius of curvature
n =  a./sqrt(1-e^2*sind(phi).^2);% prime verticl radius of curvature
%% Final Transform
x_ecef = (n + h).*cosd(phi).*cosd(lambda);
y_ecef = (n + h).*cosd(phi).*sind(lambda);
z_ecef = (n.*(1-e^2)+h).*sind(phi);

ECEF_COORDS = [x_ecef;y_ecef;z_ecef];
if data_trans
    ECEF_COORDS = ECEF_COORDS';
end