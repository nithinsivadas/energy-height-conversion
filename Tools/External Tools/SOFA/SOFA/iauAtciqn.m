%  - - - - - - - - - -
%   i a u A t c i q n
%  - - - - - - - - - -
%
%  Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
%  star-independent astrometry parameters plus a list of light-
%  deflecting bodies.
%
%  Use of this function is appropriate when efficiency is important and
%  where many star positions are to be transformed for one date.  The
%  star-independent parameters can be obtained by calling one of the
%  functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
%
%
%  If the only light-deflecting body to be taken into account is the
%  Sun, the iauAtciq function can be used instead.  If in addition the
%  parallax and proper motions are zero, the iauAtciqz function can be
%  used.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     rc,dc         ICRS RA,Dec at J2000.0 (radians)
%     pr            RA proper motion (radians/year; Note 3)
%     pd            Dec proper motion (radians/year)
%     px            parallax (arcsec)
%     rv            radial velocity (km/s, +ve if receding)
%     astrom         star-independent astrometry parameters:
%      pmt           PM time interval (SSB, Julian years)
%      eb(3)         SSB to observer (vector, au)
%      eh(3)         Sun to observer (unit vector)
%      em            distance from Sun to observer (au)
%      v(3)          barycentric observer velocity (vector, c)
%      bm1           sqrt(1-|v|^2): reciprocal of Lorenz factor
%      bpn(3,3)      bias-precession-nutation matrix
%      along         longitude + s' (radians)
%      xpl           polar motion xp wrt local meridian (radians)
%      ypl           polar motion yp wrt local meridian (radians)
%      sphi          sine of geodetic latitude
%      cphi          cosine of geodetic latitude
%      diurab        magnitude of diurnal aberration vector
%      eral          "local" Earth rotation angle (radians)
%      refa          refraction constant A (radians)
%      refb          refraction constant B (radians)
%      n             number of bodies (Note 3)
%      b(n)          data for each of the n bodies (Notes 3,4):
%      bm            mass of the body (solar masses, Note 5)
%      dl            deflection limiter (Note 6)
%      pv(2,3)       barycentric PV of the body (au, au/day)
%
%  Returned:
%     ri,di       CIRS RA,Dec (radians)
%
%  Notes:
%  1) Star data for an epoch other than J2000.0 (for example from the
%     Hipparcos catalog, which has an epoch of J1991.25) will require a
%     preliminary call to iauPmsafe before use.
%
%  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
%
%  3) The struct b contains n entries, one for each body to be
%     considered.  If n = 0, no gravitational light deflection will be
%     applied, not even for the Sun.
%
%  4) The struct b should include an entry for the Sun as well as for
%     any planet or other body to be taken into account.  The entries
%     should be in the order in which the light passes the body.
%
%  5) In the entry in the b struct for body i, the mass parameter
%     b[i].bm can, as required, be adjusted in order to allow for such
%     effects as quadrupole field.
%
%  6) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
%     the angular separation (in radians) between star and body at
%     which limiting is applied.  As phi shrinks below the chosen
%     threshold, the deflection is artificially reduced, reaching zero
%     for phi = 0.   Example values suitable for a terrestrial
%     observer, together with masses, are as follows:
%
%        body i     b[i].bm        b[i].dl
%
%        Sun        1.0            6e-6
%        Jupiter    0.00095435     3e-9
%        Saturn     0.00028574     3e-10
%
%  7) For efficiency, validation of the contents of the b array is
%     omitted.  The supplied masses must be greater than zero, the
%     position and velocity vectors must be right, and the deflection
%     limiter greater than zero.
%
%  Called:
%     iauPmpx      proper motion and parallax
%     iauLdn       light deflection by n bodies
%     iauAb        stellar aberration
%     iauRxp       product of r-matrix and pv-vector
%     iauC2s       p-vector to spherical
%     iauAnp       normalize angle into range 0 to 2pi
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ri, di] = iauAtciqn(rc, dc, pr, pd, px, rv, astrom, n, b)

% Proper motion and parallax, giving BCRS coordinate direction.
pco = iauPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb);

% Light deflection, giving BCRS natural direction.
pnat = iauLdn(n, b, astrom.eb, pco);

% Aberration, giving GCRS proper direction.
ppr = iauAb(pnat, astrom.v, astrom.em, astrom.bm1);

% Bias-precession-nutation, giving CIRS proper direction.
PI = iauRxp(astrom.bpn, ppr);

% CIRS RA,Dec.
[w, di] = iauC2s(PI);
ri = iauAnp(w);

