%  - - - - - - - -
%   i a u A p c s
%  - - - - - - - -
%
%  For an observer whose geocentric position and velocity are known,
%  prepare star-independent astrometry parameters for transformations
%  between ICRS and GCRS.  The Earth ephemeris is supplied by the
%  caller.
%
%  The parameters produced by this function are required in the space
%  motion, parallax, light deflection and aberration parts of the
%  astrometric transformation chain.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     date1         TDB as a 2-part...
%     date2         ...Julian Date (Note 1)
%     pv(2,3)       observer's geocentric pos/vel (m, m/s)
%     ebpv(2,3)     Earth barycentric PV (au, au/day)
%     ehp(3)        Earth heliocentric P (au)
%
%  Returned:
%     astrom    star-independent astrometry parameters:
%      pmt           PM time interval (SSB, Julian years)
%      eb(3)         SSB to observer (vector, au)
%      eh(3)         Sun to observer (unit vector)
%      em            distance from Sun to observer (au)
%      v(3)          barycentric observer velocity (vector, c)
%      bm1           sqrt(1-|v|^2): reciprocal of Lorenz factor
%      bpn(3,3)      bias-precession-nutation matrix
%      along         unchanged
%      xpl           unchanged
%      ypl           unchanged
%      sphi          unchanged
%      cphi          unchanged
%      diurab        unchanged
%      eral          unchanged
%      refa          unchanged
%      refb          unchanged
%
%  Notes:
%  1) The TDB date date1+date2 is a Julian Date, apportioned in any
%     convenient way between the two arguments.  For example,
%     JD(TDB)=2450123.7 could be expressed in any of these ways, among
%     others:
%
%            date1          date2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in cases
%     where the loss of several decimal digits of resolution is
%     acceptable.  The J2000 method is best matched to the way the
%     argument is handled internally and will deliver the optimum
%     resolution.  The MJD method and the date & time methods are both
%     good compromises between resolution and convenience.  For most
%     applications of this function the choice will not be at all
%     critical.
%
%     TT can be used instead of TDB without any significant impact on
%     accuracy.
%
%  2) All the vectors are with respect to BCRS axes.
%
%  3) Providing separate arguments for (i) the observer's geocentric
%     position and velocity and (ii) the Earth ephemeris is done for
%     convenience in the geocentric, terrestrial and Earth orbit cases.
%     For deep space applications it maybe more convenient to specify
%     zero geocentric position and velocity and to supply the
%     observer's position and velocity information directly instead of
%     with respect to the Earth.  However, note the different units:
%     m and m/s for the geocentric vectors, au and au/day for the
%     heliocentric and barycentric vectors.
%
%  4) In cases where the caller does not wish to provide the Earth
%     ephemeris, the function iauApcs13 can be used instead of the
%     present function.  This computes the Earth ephemeris using the
%     SOFA function iauEpv00.
%
%  5) This is one of several functions that inserts into the astrom
%     structure star-independent parameters needed for the chain of
%     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
%
%     The various functions support different classes of observer and
%     portions of the transformation chain:
%
%          functions         observer        transformation
%
%       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
%       iauApci iauApci13    terrestrial     ICRS <-> CIRS
%       iauApco iauApco13    terrestrial     ICRS <-> observed
%       iauApcs iauApcs13    space           ICRS <-> GCRS
%       iauAper iauAper13    terrestrial     update Earth rotation
%       iauApio iauApio13    terrestrial     CIRS <-> observed
%
%     Those with names ending in "13" use contemporary SOFA models to
%     compute the various ephemerides.  The others accept ephemerides
%     supplied by the caller.
%
%     The transformation from ICRS to GCRS covers space motion,
%     parallax, light deflection, and aberration.  From GCRS to CIRS
%     comprises frame bias and precession-nutation.  From CIRS to
%     observed takes account of Earth rotation, polar motion, diurnal
%     aberration and parallax (unless subsumed into the ICRS <-> GCRS
%     transformation), and atmospheric refraction.
%
%  6) The context structure astrom produced by this function is used by
%     iauAtciq* and iauAticq*.
%
%  Called:
%     iauCp        copy p-vector
%     iauPm        modulus of p-vector
%     iauPn        decompose p-vector into modulus and direction
%     iauIr        initialize r-matrix to identity
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function astrom = iauApcs(date1, date2, pv, ebpv, ehp, astrom)

constants

% au/d to m/s
AUDMS = DAU/DAYSEC;

% Light time for 1 AU (day)
CR = AULT/DAYSEC;

pb = zeros(3,1);
vb = zeros(3,1);
ph = zeros(3,1);

% Time since reference epoch, years (for proper motion calculation).
astrom.pmt = ( (date1 - DJ00) + date2 ) / DJY;

% Adjust Earth ephemeris to observer.
for i = 1:3
      dp = pv(1,i) / DAU;
      dv = pv(2,i) / AUDMS;
      pb(i) = ebpv(1,i) + dp;
      vb(i) = ebpv(2,i) + dv;
      ph(i) = ehp(i) + dp;
end

% Barycentric position of observer (au).
astrom.eb = iauCp(pb);

% Heliocentric direction and distance (unit vector and au).
[astrom.em, astrom.eh] = iauPn(ph);

% Barycentric vel. in units of c, and reciprocal of Lorenz factor.
v2 = 0;
for i = 1:3
      w = vb(i) * CR;
      astrom.v(i) = w;
      v2 = v2 + w*w;
end
astrom.bm1 = sqrt(1 - v2);

% Reset the NPB matrix.
astrom.bpn = zeros(3);
astrom.bpn = iauIr(astrom.bpn);

