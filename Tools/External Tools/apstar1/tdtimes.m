function [ttjd, secdif] = tdtimes (tdbjd)

% this function computes the terrestrial time (tt)
% julian date corresponding to a barycentric dynamical
% time (tdb) julian date. expressions used in this
% version are approximations resulting in accuracies
% of about 20 microseconds. see explanatory supplement
% to the astronomical almanac, pp. 42-44 and 316.

% input

%  tdbjd = tdb julian date

% output

%  ttjd = tt julian date

%  secdif = difference tdbjd-ttjd, in seconds (out)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seccon = 206264.8062470964;

rev = 1296000;

% tdb julian date of epoch j2000.0

t0 = 2451545;

% eccentricity of earth-moon barycenter orbit

ecc = 0.01671022;

tdays = tdbjd - t0;

m = (357.51716 + 0.985599987 * tdays) * 3600;

l = (280.46435 + 0.985609100 * tdays) * 3600;

lj = (34.40438 + 0.083086762 * tdays) * 3600;

m  = mod(m, rev) / seccon;

l  = mod (l, rev) / seccon;

lj = mod (lj, rev) / seccon;

e  = m + ecc * sin(m) + 0.5 * ecc^2 * sin (2 * m);

secdif = 1.658e-3 * sin(e) + 20.73e-6 * sin(l - lj);

ttjd = tdbjd - secdif / 86400;

