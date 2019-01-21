function [pos, vel] = terra (glon, glat, ht, st)

% this functions computes the position and velocity vectors of
% a terrestrial observer with respect to the center of the earth.

% input

%  glon = longitude of observer with respect to reference
%         meridian (east +) in degrees
%  glat = geodetic latitude (north +) of observer in degrees
%  ht   = height of observer in meters
%  st   = local apparent sidereal time at reference meridian
%         in hours

% output

%  pos = position vector of observer with respect to center
%        of earth, equatorial rectangular coordinates,
%        referred to true equator and equinox of date,
%        components in au
%  vel = velocity vector of observer with respect to center
%        of earth, equatorial rectangular coordinates,
%        referred to true equator and equinox of date,
%        components in au/day

%  note: if reference meridian is greenwich and st=0, pos
%        is effectively referred to equator and greenwich

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global stlocl

seccon = 206264.8062470964;

erad = 6378.140;

f = 0.00335281;

omega = 7.2921151467e-5;

kmau = 1.49597870e8;

% compute parameters relating to geodetic to geocentric conversion

df2 = (1 - f)^2;

phi = glat * 3600 / seccon;

sinphi = sin(phi);

cosphi = cos(phi);

c = 1 / sqrt(cosphi^2 + df2 * sinphi^2);

s = df2 * c;

ach = erad * c + ht / 1000;

ash = erad * s + ht / 1000;

% compute local sidereal time factors

stlocl = (st * 54000 + glon * 3600) / seccon;

sinst = sin(stlocl);

cosst = cos(stlocl);

% compute position vector components in km

pos(1) = ach * cosphi * cosst;

pos(2) = ach * cosphi * sinst;

pos(3) = ash * sinphi;

% compute velocity vector components in km/sec

vel(1) = -omega * ach * cosphi * sinst;

vel(2) =  omega * ach * cosphi * cosst;

vel(3) =  0;

% convert position and velocity components to au and au/day

for j = 1:1:3
    pos(j) = pos(j) / kmau;
    vel(j) = vel(j) / kmau * 86400;
end

