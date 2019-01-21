%  - - - - - - - - - -
%   i a u A t c i q z
%  - - - - - - - - - -
%
%  Quick ICRS to CIRS transformation, given precomputed star-
%  independent astrometry parameters, and assuming zero parallax and
%  proper motion.
%
%  Use of this function is appropriate when efficiency is important and
%  where many star positions are to be transformed for one date.  The
%  star-independent parameters can be obtained by calling one of the
%  functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
%
%  The corresponding function for the case of non-zero parallax and
%  proper motion is iauAtciq.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     rc,dc       ICRS astrometric RA,Dec (radians)
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
%
%  Returned:
%     ri,di       CIRS RA,Dec (radians)
%
%  Note:
%     All the vectors are with respect to BCRS axes.
%
%  References:
%     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
%     the Astronomical Almanac, 3rd ed., University Science Books
%     (2013).
%     Klioner, Sergei A., "A practical relativistic model for micro-
%     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
%
%  Called:
%     iauS2c       spherical coordinates to unit vector
%     iauLdsun     light deflection due to Sun
%     iauAb        stellar aberration
%     iauRxp       product of r-matrix and p-vector
%     iauC2s       p-vector to spherical
%     iauAnp       normalize angle into range +/- pi
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ri, di] = iauAtciqz(rc, dc, astrom)

% BCRS coordinate direction (unit vector).
pco = iauS2c(rc, dc);

% Light deflection by the Sun, giving BCRS natural direction.
pnat = iauLdsun(pco, astrom.eh, astrom.em);

% Aberration, giving GCRS proper direction.
ppr = iauAb(pnat, astrom.v, astrom.em, astrom.bm1);

% Bias-precession-nutation, giving CIRS proper direction.
PI = iauRxp(astrom.bpn, ppr);

% CIRS RA,Dec.
[w, di] = iauC2s(PI);
ri = iauAnp(w);

