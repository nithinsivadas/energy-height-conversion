%  - - - - - - - -
%   i a u A p e r
%  - - - - - - - -
%
%  In the star-independent astrometry parameters, update only the
%  Earth rotation angle, supplied by the caller explicitly.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     theta          Earth rotation angle (radians, Note 2)
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
%  1) This function exists to enable sidereal-tracking applications to
%     avoid wasteful recomputation of the bulk of the astrometry
%     parameters:  only the Earth rotation is updated.
%
%  2) For targets expressed as equinox based positions, such as
%     classical geocentric apparent (RA,Dec), the supplied theta can be
%     Greenwich apparent sidereal time rather than Earth rotation
%     angle.
%
%  3) The function iauAper13 can be used instead of the present
%     function, and starts from UT1 rather than ERA itself.
%
%  4) This is one of several functions that inserts into the astrom
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
%  This revision:   2013 September 25
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function astrom = iauAper(theta, astrom)

astrom.eral = theta + astrom.along;

