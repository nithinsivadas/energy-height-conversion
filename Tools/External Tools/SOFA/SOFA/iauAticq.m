%  - - - - - - - - -
%   i a u A t i c q
%  - - - - - - - - -
%
%  Quick CIRS RA,Dec to ICRS astrometric place, given the star-
%  independent astrometry parameters.
%
%  Use of this function is appropriate when efficiency is important and
%  where many star positions are all to be transformed for one date.
%  The star-independent astrometry parameters can be obtained by
%  calling one of the functions iauApci[13], iauApcg[13], iauApco[13]
%  or iauApcs[13].
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     ri,di          CIRS RA,Dec (radians)
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
%     rc,dc       ICRS astrometric RA,Dec (radians)
%
%  Notes:
%
%  1) Only the Sun is taken into account in the light deflection
%     correction.
%
%  2) Iterative techniques are used for the aberration and light
%     deflection corrections so that the functions iauAtic13 (or
%     iauAticq) and iauAtci13 (or iauAtciq) are accurate inverses;
%     even at the edge of the Sun's disk the discrepancy is only about
%     1 nanoarcsecond.
%
%  Called:
%     iauS2c       spherical coordinates to unit vector
%     iauTrxp      product of transpose of r-matrix and p-vector
%     iauZp        zero p-vector
%     iauAb        stellar aberration
%     iauLdsun     light deflection by the Sun
%     iauC2s       p-vector to spherical
%     iauAnp       normalize angle into range +/- pi
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rc, dc] = iauAticq(ri, di, astrom)

before = zeros(3,1);
pnat = zeros(3,1);
pco = zeros(3,1);

% CIRS RA,Dec to Cartesian.
PI = iauS2c(ri, di);

% Bias-precession-nutation, giving GCRS proper direction.
ppr = iauTrxp(astrom.bpn, PI);

% Aberration, giving GCRS natural direction.
d = zeros(3,1);
d = iauZp(d);
for j = 1:2
    r2 = 0;
    for i = 1:3
        w = ppr(i) - d(i);
        before(i) = w;
        r2 = r2 + w*w;
    end
    r = sqrt(r2);
    for i = 1:3
        before(i) = before(i) / r;
    end
    after = iauAb(before, astrom.v, astrom.em, astrom.bm1);
    r2 = 0;
    for i = 1:3
        d(i) = after(i) - before(i);
        w = ppr(i) - d(i);
        pnat(i) = w;
        r2 = r2 + w*w;
    end
    r = sqrt(r2);
    for i = 1:3
        pnat(i) = pnat(i) / r;
    end
end

% Light deflection by the Sun, giving BCRS coordinate direction.
iauZp(d);
for j = 1:5
      r2 = 0;
      for i = 1:3
         w = pnat(i) - d(i);
         before(i) = w;
         r2 = r2 + w*w;
      end
      r = sqrt(r2);
      for i = 1:3
         before(i) = before(i) / r;
      end
      after = iauLdsun(before, astrom.eh, astrom.em);
      r2 = 0;
      for i = 1:3
         d(i) = after(i) - before(i);
         w = pnat(i) - d(i);
         pco(i) = w;
         r2 = r2 + w*w;
      end
      r = sqrt(r2);
      for i = 1:3
         pco(i) = pco(i) / r;
      end
end

% ICRS astrometric RA,Dec.
[w, dc] = iauC2s(pco);
rc = iauAnp(w);

