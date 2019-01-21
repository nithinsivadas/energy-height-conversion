%  - - - - - - - - -
%   i a u A t c i q
%  - - - - - - - - -
%
%  Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
%  star-independent astrometry parameters.
%
%  Use of this function is appropriate when efficiency is important and
%  where many star positions are to be transformed for one date.  The
%  star-independent parameters can be obtained by calling one of the
%  functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
%
%  If the parallax and proper motions are zero the iauAtciqz function
%  can be used instead.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     rc,dc       ICRS RA,Dec at J2000.0 (radians)
%     pr          RA proper motion (radians/year; Note 3)
%     pd          Dec proper motion (radians/year)
%     px          parallax (arcsec)
%     rv          radial velocity (km/s, +ve if receding)
%     astrom      star-independent astrometry parameters:
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
%
%  Returned:
%     ri,di       CIRS RA,Dec (radians)
%
%  Notes:
%  1) All the vectors are with respect to BCRS axes.
%
%  2) Star data for an epoch other than J2000.0 (for example from the
%     Hipparcos catalog, which has an epoch of J1991.25) will require a
%     preliminary call to iauPmsafe before use.
%
%  3) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
%
%  Called:
%     iauPmpx      proper motion and parallax
%     iauLdsun     light deflection by the Sun
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
function [ri, di] = iauAtciq(rc, dc, pr, pd, px, rv, astrom)

% Proper motion and parallax, giving BCRS coordinate direction.
pco = iauPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb);

% Light deflection by the Sun, giving BCRS natural direction.
pnat = iauLdsun(pco, astrom.eh, astrom.em);

% Aberration, giving GCRS proper direction.
ppr = iauAb(pnat, astrom.v, astrom.em, astrom.bm1);

% Bias-precession-nutation, giving CIRS proper direction.
PI = iauRxp(astrom.bpn, ppr);

% CIRS RA,Dec.
[w, di] = iauC2s(PI);
ri = iauAnp(w);

