function WGS_COORDS = ecef2wgs(ECEF_COORDS)
% ecef2wgs.m
% ECEF_COORDS = wgs2ecef(WGS_COORDS)
% by John Swoboda 1/8/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
% This function takes a set of X,Y Z ECEF coordinates and puts them into
% WGS84 coordinate of latitude longitude and height from the surface of the 
% earth.  The input is assumed to be a 3xN (or Nx3) matrix of ECEF 
% coordinates in meters.  It then outputs a 3xN matrix of WGS coordinates 
% in degrees and meters.  This has been compared to matlab's internal
% ecef2geodetic and has comparible results with a normalized difference on 
% the order of 10^-16. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% ECEF_COORDS - A 3xN matrix (can also be a Nx3 matrix) with X, Y and Z in   
% the ECEF coordinate space.  Coordinates are in units of meters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output  
% WGS_COORDS - A 3xN matrix  with the latitude longitude and height 
% coordinates.  The matrix is broken up in the following way, the first row 
% is latitude in degrees, the second row is longitude in degrees and the 
% last row is height from the surface of the earth in meters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference
% Vermeille H (2002) Direct transformation from geocentric coordinates
% to geodetic coordinates, Journal of Geodesy, vol. 76, no. 8, pp 451 - 454,
% Nov 2002, http://link.springer.com/article/10.1007%2Fs00190-002-0273-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input
[n,m] = size(ECEF_COORDS);
if n~=3&&m==3
    WGS_COORDS = ECEF_COORDS.';
elseif n==3&&m==3
    warning('Both dimensions have 3 elements, assuming lat, long and height are along rows');
elseif n~=3 &&m~=3
    error('Neither dimension has a size of 3');
end
%% Get Lat, Long and Height
X = ECEF_COORDS(1,:);
Y = ECEF_COORDS(2,:);
Z = ECEF_COORDS(3,:);
%% Set the constants
a = 6378137; % semi-major axis in meters
f = 1/298.257223563; % the flattening factor
b = a*(1-f);% semiminor axis in meters 
e = sqrt(a^2-b^2)/a;% first eccentricity

%% Algorithm
% this is taken straight from the Vermeille paper
p = (X.^2+Y.^2)/a^2;

q = ((1-e^2)/a^2)*Z.^2;

r = (p+q-e^4)/6;

s = e^4*(p.*q)./(4.*r.^3);

t = nthroot(1+s+sqrt(s.*(2+s)),3);

u = r.*(1+t+(1./t));

v = sqrt(u.^2+(e^4)*q);

w = e^2*((u+v-q)./(2*v));

k = sqrt(u+v+w.^2)-w;

D = k.*sqrt(X.^2+Y.^2)./(k+e^2);

%% Final Form
% use the atan2 function for more numerical stability
long = atan2(Y,X)*180/pi;

lat = atan2(Z,D)*180/pi;

h = ((k+e^2-1)./(k)).*sqrt(D.^2+Z.^2);
% put into the final form
WGS_COORDS = [lat;long;h];