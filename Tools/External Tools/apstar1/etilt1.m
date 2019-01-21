function [oblm, oblt, eqeq, dpsi, deps] = etilt1 (tjd)

% this function computes quantities related to the orientation
% of the earth's rotation axis at julian date tjd

% nutation parameters from a jpl binary ephemeris

% input

%  tjd = tdb julian date for orientation parameters

% output

%  oblm = mean obliquity of the ecliptic at date tjd (degrees)
%  oblt = true obliquity of the ecliptic at date tjd (degrees)
%  eqeq = equation of the equinoxes at date tjd (arc seconds)
%  dpsi = nutation in longitude at date tjd (arc seconds)
%  deps = nutation in obliquity at date tjd (arc seconds)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global tlaste obm obt ee psi eps

t0 = 2451545;

seccon = 206264.8062470964;

if (abs(tjd - tlaste) < 1e-6)
   oblm = obm;
   oblt = obt;
   eqeq = ee;
   dpsi = psi;
   deps = eps;
   return;
end

t = (tjd - t0) / 36525;

t2 = t * t;

t3 = t2 * t;

% obtain nutation parameters in seconds of arc

rrd = jplephem(tjd, 14, 0);

psi = rrd(1) * seccon;

eps = rrd(2) * seccon;

% compute mean obliquity of the ecliptic in seconds of arc

obm = 84381.4480 - 46.8150 * t - 0.00059 * t2 + 0.001813 * t3;

% compute true obliquity of the ecliptic in seconds of arc

obt = obm + eps;

% compute equation of the equinoxes in seconds of time

ee = psi / 15 * cos (obt/seccon);

% convert obliquity values to degrees

obm = obm / 3600;

obt = obt / 3600;

tlaste = tjd;

oblm = obm;
oblt = obt;
eqeq = ee;
dpsi = psi;
deps = eps;
