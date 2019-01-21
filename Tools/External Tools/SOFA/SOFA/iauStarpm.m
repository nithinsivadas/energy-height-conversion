%  - - - - - - - - - -
%   i a u S t a r p m
%  - - - - - - - - - -
%
%  Star proper motion:  update star catalog data for space motion.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     ra1        right ascension (radians), before
%     dec1       declination (radians), before
%     pmr1       RA proper motion (radians/year), before
%     pmd1       Dec proper motion (radians/year), before
%     px1        parallax (arcseconds), before
%     rv1        radial velocity (km/s, +ve = receding), before
%     ep1a       "before" epoch, part A (Note 1)
%     ep1b       "before" epoch, part B (Note 1)
%     ep2a       "after" epoch, part A (Note 1)
%     ep2b       "after" epoch, part B (Note 1)
%
%  Returned:
%     ra2        right ascension (radians), after
%     dec2       declination (radians), after
%     pmr2       RA proper motion (radians/year), after
%     pmd2       Dec proper motion (radians/year), after
%     px2        parallax (arcseconds), after
%     rv2        radial velocity (km/s, +ve = receding), after
%
%  Notes:
%  1) The starting and ending TDB dates ep1a+ep1b and ep2a+ep2b are
%     Julian Dates, apportioned in any convenient way between the two
%     parts (A and B).  For example, JD(TDB)=2450123.7 could be
%     expressed in any of these ways, among others:
%
%             epna          epnb
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable.  The J2000 method is best matched to the way
%     the argument is handled internally and will deliver the
%     optimum resolution.  The MJD method and the date & time methods
%     are both good compromises between resolution and convenience.
%
%  2) In accordance with normal star-catalog conventions, the object's
%     right ascension and declination are freed from the effects of
%     secular aberration.  The frame, which is aligned to the catalog
%     equator and equinox, is Lorentzian and centered on the SSB.
%
%     The proper motions are the rate of change of the right ascension
%     and declination at the catalog epoch and are in radians per TDB
%     Julian year.
%
%     The parallax and radial velocity are in the same frame.
%
%  3) Care is needed with units.  The star coordinates are in radians
%     and the proper motions in radians per Julian year, but the
%     parallax is in arcseconds.
%
%  4) The RA proper motion is in terms of coordinate angle, not true
%     angle.  If the catalog uses arcseconds for both RA and Dec proper
%     motions, the RA proper motion will need to be divided by cos(Dec)
%     before use.
%
%  5) Straight-line motion at constant speed, in the inertial frame,
%     is assumed.
%
%  6) An extremely small (or zero or negative) parallax is interpreted
%     to mean that the object is on the "celestial sphere", the radius
%     of which is an arbitrary (large) value (see the iauStarpv
%     function for the value used).  When the distance is overridden in
%     this way, the status, initially zero, has 1 added to it.
%
%  7) If the space velocity is a significant fraction of c (see the
%     constant VMAX in the function iauStarpv),  it is arbitrarily set
%     to zero.  When this action occurs, 2 is added to the status.
%
%  8) The relativistic adjustment carried out in the iauStarpv function
%     involves an iterative calculation.  If the process fails to
%     converge within a set number of iterations, 4 is added to the
%     status.
%
%  Called:
%     iauStarpv    star catalog data to space motion pv-vector
%     iauPvu       update a pv-vector
%     iauPdp       scalar product of two p-vectors
%     iauPvstar    space motion pv-vector to star catalog data
%
%  This revision:  2008 May 16
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ra2, dec2, pmr2, pmd2, px2, rv2] = iauStarpm(ra1, dec1, pmr1, ...
          pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b)

constants

% RA,Dec etc. at the "before" epoch to space motion pv-vector.
pv1 = iauStarpv(ra1, dec1, pmr1, pmd1, px1, rv1);

% Light time when observed (days).
tl1 = iauPm(pv1(1,:)) / DC;

% Time interval, "before" to "after" (days).
dt = (ep2a - ep1a) + (ep2b - ep1b);

% Move star along track from the "before" observed position to the
% "after" geometric position.
pv = iauPvu(dt + tl1, pv1);

% From this geometric position, deduce the observed light time (days)
% at the "after" epoch (with theoretically unneccessary error check).
r2 = iauPdp(pv(1,:), pv(1,:));
rdv = iauPdp(pv(1,:), pv(2,:));
v2 = iauPdp(pv(2,:), pv(2,:));
c2mv2 = DC*DC - v2;
tl2 = (-rdv + sqrt(rdv*rdv + c2mv2*r2)) / c2mv2;

% Move the position along track from the observed place at the
% "before" epoch to the observed place at the "after" epoch.
pv2 = iauPvu(dt + (tl1 - tl2), pv1);

% Space motion pv-vector to RA,Dec etc. at the "after" epoch.
[ra2, dec2, pmr2, pmd2, px2, rv2] = iauPvstar(pv2);

