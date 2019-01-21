%  - - - - - - - -
%   i a u A p i o
%  - - - - - - - -
%
%  For a terrestrial observer, prepare star-independent astrometry
%  parameters for transformations between CIRS and observed
%  coordinates.  The caller supplies the Earth orientation information
%  and the refraction constants as well as the site coordinates.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     sp           the TIO locator s' (radians, Note 1)
%     theta        Earth rotation angle (radians)
%     elong        longitude (radians, east +ve, Note 2)
%     phi          geodetic latitude (radians, Note 2)
%     hm           height above ellipsoid (m, geodetic Note 2)
%     xp,yp        polar motion coordinates (radians, Note 3)
%     refa         refraction constant A (radians, Note 4)
%     refb         refraction constant B (radians, Note 4)
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
%      along         longitude + s' (radians)
%      xpl           polar motion xp wrt local meridian (radians)
%      ypl           polar motion yp wrt local meridian (radians)
%      sphi          sine of geodetic latitude
%      cphi          cosine of geodetic latitude
%      diurab        magnitude of diurnal aberration vector
%      eral          "local" Earth rotation angle (radians)
%      refa          refraction constant A (radians)
%      refb          refraction constant B (radians)
%
%  Notes:
%  1) sp, the TIO locator s', is a tiny quantity needed only by the
%     most precise applications.  It can either be set to zero or
%     predicted using the SOFA function iauSp00.
%
%  2) The geographical coordinates are with respect to the WGS84
%     reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
%     longitude required by the present function is east-positive
%     (i.e. right-handed), in accordance with geographical convention.
%
%  3) The polar motion xp,yp can be obtained from IERS bulletins.  The
%     values are the coordinates (in radians) of the Celestial
%     Intermediate Pole with respect to the International Terrestrial
%     Reference System (see IERS Conventions 2003), measured along the
%     meridians 0 and 90 deg west respectively.  For many applications,
%     xp and yp can be set to zero.
%
%     Internally, the polar motion is stored in a form rotated onto the
%     local meridian.
%
%  4) The refraction constants refa and refb are for use in a
%     dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
%     (i.e. refracted) zenith distance and dZ is the amount of
%     refraction.
%
%  5) It is advisable to take great care with units, as even unlikely
%     values of the input parameters are accepted and processed in
%     accordance with the models used.
%
%  6) In cases where the caller does not wish to provide the Earth
%     rotation information and refraction constants, the function
%     iauApio13 can be used instead of the present function.  This
%     starts from UTC and weather readings etc. and computes suitable
%     values using other SOFA functions.
%
%  7) This is one of several functions that inserts into the astrom
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
%  8) The context structure astrom produced by this function is used by
%     iauAtioq and iauAtoiq.
%
%  Called:
%     iauPvtob     position/velocity of terrestrial station
%     iauAper      astrometry parameters: update ERA
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function astrom = iauApio(sp, theta, elong, phi, hm, xp, yp, refa, refb, astrom)

constants

% Longitude with adjustment for TIO locator s'.
astrom.along = elong + sp;

% Polar motion, rotated onto the local meridian.
sl = sin(astrom.along);
cl = cos(astrom.along);
astrom.xpl = xp*cl - yp*sl;
astrom.ypl = xp*sl + yp*cl;

% Functions of latitude.
astrom.sphi = sin(phi);
astrom.cphi = cos(phi);

% Observer's geocentric position and velocity (m, m/s, CIRS).
pv = iauPvtob(elong, phi, hm, xp, yp, sp, theta);

% Magnitude of diurnal aberration vector.
astrom.diurab = sqrt(pv(2,1)*pv(2,1)+pv(2,2)*pv(2,2)) / CMPS;

% Refraction constants.
astrom.refa = refa;
astrom.refb = refb;

% Local Earth rotation angle.
astrom = iauAper(theta, astrom);

