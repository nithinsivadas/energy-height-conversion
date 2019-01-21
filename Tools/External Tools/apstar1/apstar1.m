function [ra, dec] = apstar1 (tjd, ujd, n, topo, glon, glat, ht, ...
                              ram, decm, pmra, pmdec, parlax, radvel)

% this subroutine computes the geocentric or topocentric apparent place
% of a star, given its mean place, proper motion, parallax, and radial
% velocity for j2000.0. the coordinates of the earth are determined
% using jpleph.

% input

%  tjd    = tdt julian date for apparent place
%  ujd    = ut1 julian date for apparent topocentric place
%  n      = body identification number for the earth
%  topo   = type of apparent place calculation
%         = 0 ==> geocentric
%         = 1 ==> topocentric
%  ram    = mean right ascension j2000.0 (hours)
%  decm   = mean declination j2000.0 (degrees)
%  pmra   = proper motion in ra (seconds of time/julian century)
%  pmdec  = proper motion in dec (seconds of arc/julian century)
%  parlax = parallax (seconds of arc)
%  radvel = radial velocity (kilometers per second)

% output

%  ra  = apparent geocentric or topocentric right ascension,
%        referred to true equator and equinox of date (hours)
%  dec = apparent geocentric or topocentric declination,
%        referred to true equator and equinox of date (degrees)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 173.1446333;

t0 = 2451545;

% compute t1, the tdb julian date corresponding to tjd

[x, secdif] = tdtimes (tjd);

t1 = tjd + secdif / 86400;

% get position and velocity of the earth wrt barycenter
% of solar system and wrt center of sun

[pebs, vebs, ierr] = solsys (t1, n, 0);

if (ierr ~= 0)
   ra = 0;
   dec = 0;
   return;
end

[pess, vess, ierr] = solsys (t1, n, 1);

if (ierr ~= 0)
   ra = 0;
   dec = 0;
   return;
end

for j = 1:1:3
    pb(j) = pebs(j);
    vb(j) = vebs(j);
    ps(j) = pess(j);
    vs(j) = vess(j);
end

rm = ram;

dm = decm;

pmr = pmra;

pmd = pmdec;

pi = parlax;

rv = radvel;

if (topo == 1)
   % get position and velocity of observer wrt center of earth

   st = gast2 (ujd, 0, 0);

   [x1, x2, eqeq, x3, x4] = etilt1 (t1);

   gast = st + eqeq / 3600;

   [pos1, vel1] = terra (glon, glat, ht, gast);

   pos2 = nutate1 (-t1, pos1);

   pog = precess (t1, pos2, t0);

   vel2 = nutate1 (-t1, vel1);

   vog = precess (t1, vel2, t0);

   % compute position and velocity of observer wrt barycenter of
   % solar system and wrt center of sun

   for j = 1:1:3
       pb(j) = pebs(j) + pog(j);
       vb(j) = vebs(j) + vog(j);
       ps(j) = pess(j) + pog(j);
       vs(j) = vess(j) + vog(j);
    end
end
 
[pos1, vel1] = vectrs(rm, dm, pmr, pmd, pi, rv);

pos2 = propmo(t0, pos1, vel1, t1);

[pos3, tlight] = geocen (pos2, pb);

pos4 = sunfld (pos3, ps);

pos5 = aberat (pos4, vb, tlight);

pos6 = precess (t0, pos5, t1);

pos7 = nutate1 (t1, pos6);

[r, d] = angles (pos6);

ra = r;

dec = d;


