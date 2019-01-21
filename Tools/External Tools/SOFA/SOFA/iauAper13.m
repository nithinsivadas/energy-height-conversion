%  - - - - - - - - - -
%   i a u A p e r 1 3
%  - - - - - - - - - -
%
%  In the star-independent astrometry parameters, update only the
%  Earth rotation angle.  The caller provides UT1, (n.b. not UTC).
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     ut11          UT1 as a 2-part...
%     ut12          ...Julian Date (Note 1)
%     astrom         star-independent astrometry parameters:
%      pmt           not used
%      eb(3)         not used
%      eh(3)         not used
%      em            not used
%      v(3)          not used
%      bm1           not used
%      bpn(3,3)      not used
%      along         longitude + s' (radians)
%      xpl           not used
%      ypl           not used
%      sphi          not used
%      cphi          not used
%      diurab        not used
%      eral          not used
%      refa          not used
%      refb          not used
%
%  Returned:
%     astrom         star-independent astrometry parameters:
%      pmt           unchanged
%      eb(3)         unchanged
%      eh(3)         unchanged
%      em            unchanged
%      v(3)          unchanged
%      bm1           unchanged
%      bpn(3,3)      unchanged
%      along         unchanged
%      xpl           unchanged
%      ypl           unchanged
%      sphi          unchanged
%      cphi          unchanged
%      diurab        unchanged
%      eral          "local" Earth rotation angle (radians)
%      refa          unchanged
%      refb          unchanged
%
%  Notes:
%  1) The UT1 date (n.b. not UTC) ut11+ut12 is a Julian Date,
%     apportioned in any convenient way between the arguments ut11 and
%     ut12.  For example, JD(UT1)=2450123.7 could be expressed in any
%     of these ways, among others:
%
%            ut11           ut12
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in cases
%     where the loss of several decimal digits of resolution is
%     acceptable.  The J2000 and MJD methods are good compromises
%     between resolution and convenience.  The date & time method is
%     best matched to the algorithm used:  maximum precision is
%     delivered when the ut11 argument is for 0hrs UT1 on the day in
%     question and the ut12 argument lies in the range 0 to 1, or vice
%     versa.
%
%  2) If the caller wishes to provide the Earth rotation angle itself,
%     the function iauAper can be used instead.  One use of this
%     technique is to substitute Greenwich apparent sidereal time and
%     thereby to support equinox based transformations directly.
%
%  3) This is one of several functions that inserts into the astrom
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
%  Called:
%     iauAper      astrometry parameters: update ERA
%     iauEra00     Earth rotation angle, IAU 2000
%
%  This revision:   2013 September 25
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function astrom = iauAper13(ut11, ut12, astrom)

astrom = iauAper(iauEra00(ut11,ut12), astrom);

