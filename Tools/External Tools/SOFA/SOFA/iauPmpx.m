%  - - - - - - - -
%   i a u P m p x
%  - - - - - - - -
%
%  Proper motion and parallax.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     rc,dc       ICRS RA,Dec at catalog epoch (radians)
%     pr          RA proper motion (radians/year; Note 1)
%     pd          Dec proper motion (radians/year)
%     px          parallax (arcsec)
%     rv          radial velocity (km/s, +ve if receding)
%     pmt         proper motion time interval (SSB, Julian years)
%     pob(3)      SSB to observer vector (au)
%
%  Returned:
%     pco(3)      coordinate direction (BCRS unit vector)
%
%  Notes:
%  1) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
%
%  2) The proper motion time interval is for when the starlight
%     reaches the solar system barycenter.
%
%  3) To avoid the need for iteration, the Roemer effect (i.e. the
%     small annual modulation of the proper motion coming from the
%     changing light time) is applied approximately, using the
%     direction of the star at the catalog epoch.
%
%  References:
%     1984 Astronomical Almanac, pp B39-B41.
%     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
%     the Astronomical Almanac, 3rd ed., University Science Books
%     (2013), Section 7.2.
%
%  Called:
%     iauPdp       scalar product of two p-vectors
%     iauPn        decompose p-vector into modulus and direction
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pco = iauPmpx(rc, dc, pr, pd, px, rv, pmt, pob)

constants

% Km/s to au/year
VF = DAYSEC*DJM/DAU;

% Light time for 1 au, Julian years
AULTY = AULT/DAYSEC/DJY;

% int i;
% sr, cr, sd, cd, x, y, z, p(3), dt, pxr, w, pdz, pm(3);

% Spherical coordinates to unit vector (and useful functions).
sr = sin(rc);
cr = cos(rc);
sd = sin(dc);
cd = cos(dc);
p  = zeros(3,1);
x = cr*cd;
p(1) = x;
y = sr*cd;
p(2) = y;
z = sd;
p(3) = z;

% Proper motion time interval (y) including Roemer effect.
dt = pmt + iauPdp(p,pob)*AULTY;

% Space motion (radians per year).
pxr = px * DAS2R;
w = VF * rv * pxr;
pdz = pd * z;
pm(1) = - pr*y - pdz*cr + w*x;
pm(2) =   pr*x - pdz*sr + w*y;
pm(3) =   pd*cd + w*z;

% Coordinate direction of star (unit vector, BCRS).
for i = 1:3
    p(i) = p(i) + dt*pm(i) - pxr*pob(i);
end
[w, pco] = iauPn(p);

