function [pos, vel] = vectrs (ra, dec, pmra, pmdec, parllx, rv)

% this function converts angular quantities to vectors

% input

%  ra     = right ascension (hours)
%  dec    = declination (degrees)
%  pmra   = proper motion in ra (seconds of time per julian century)
%  pmdec  = proper motion in dec per julian century (seconds of arc)
%  parllx = parallax (seconds of arc)
%  rv     = radial velocity (kilometers per second)

% output

%  pos = position vector, equatorial rectangular coordinates (au)
%  vel = velocity vector, equatorial rectangular coordinates (au/day)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seccon = 206264.8062470964;

kmau = 1.49597870e8;

% if parallax is unknown, undetermined, or zero, set it to 1e-7
% second of arc, corresponding to a distance of 10 megaparsecs

paralx = parllx;

if (paralx <= 0) 
   paralx = 1e-7;
end

% convert right ascension, declination, and parallax to
% position vector in equatorial system with units of au

dist = seccon / paralx;

r = ra * 54000 / seccon;

d = dec * 3600 / seccon;

cra = cos(r);
sra = sin(r);
cdc = cos(d);
sdc = sin(d);

pos(1) = dist * cdc * cra;
pos(2) = dist * cdc * sra;
pos(3) = dist * sdc;

% convert proper motion and radial velocity to orthogonal
% components of motion with units of au/day

pmr = pmra * 150 * cdc / (paralx * 36525);

pmd = pmdec / (paralx * 36525);

rvl = rv * 86400 / kmau;

% transform motion vector to equatorial system

vel(1) = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra;
vel(2) = pmr * cra - pmd * sdc * sra + rvl * cdc * sra;
vel(3) = pmd * cdc + rvl * sdc;

