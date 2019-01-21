%  - - - - - - - - -
%   i a u A t i c q n
%  - - - - - - - - -
%
%  Quick CIRS to ICRS astrometric place transformation, given the star-
%  independent astrometry parameters plus a list of light-deflecting
%  bodies.
%
%  Use of this function is appropriate when efficiency is important and
%  where many star positions are all to be transformed for one date.
%  The star-independent astrometry parameters can be obtained by
%  calling one of the functions iauApci[13], iauApcg[13], iauApco[13]
%  or iauApcs[13].
%
%  If the only light-deflecting body to be taken into account is the
%  Sun, the iauAticq function can be used instead.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     ri,di        CIRS RA,Dec (radians)
%     astrom       star-independent astrometry parameters:
%      pmt         PM time interval (SSB, Julian years)
%      eb(3)       SSB to observer (vector, au)
%      eh(3)       Sun to observer (unit vector)
%      em          distance from Sun to observer (au)
%      v(3)        barycentric observer velocity (vector, c)
%      bm1         sqrt(1-|v|^2): reciprocal of Lorenz factor
%      bpn(3,3)    bias-precession-nutation matrix
%      along       longitude + s' (radians)
%      xpl         polar motion xp wrt local meridian (radians)
%      ypl         polar motion yp wrt local meridian (radians)
%      sphi        sine of geodetic latitude
%      cphi        cosine of geodetic latitude
%      diurab      magnitude of diurnal aberration vector
%      eral        "local" Earth rotation angle (radians)
%      refa        refraction constant A (radians)
%      refb        refraction constant B (radians)
%      n           number of bodies (Note 3)
%      b           data for each of the n bodies (Notes 3,4):
%      bm          mass of the body (solar masses, Note 5)
%      dl          deflection limiter (Note 6)
%      pv(2,3)     barycentric PV of the body (au, au/day)
%
%  Returned:
%     rc,dc       ICRS astrometric RA,Dec (radians)
%
%  Notes:
%  1) Iterative techniques are used for the aberration and light
%     deflection corrections so that the functions iauAticqn and
%     iauAtciqn are accurate inverses; even at the edge of the Sun's
%     disk the discrepancy is only about 1 nanoarcsecond.
%
%  2) If the only light-deflecting body to be taken into account is the
%     Sun, the iauAticq function can be used instead.
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
%     iauS2c       spherical coordinates to unit vector
%     iauTrxp      product of transpose of r-matrix and p-vector
%     iauZp        zero p-vector
%     iauAb        stellar aberration
%     iauLdn       light deflection by n bodies
%     iauC2s       p-vector to spherical
%     iauAnp       normalize angle into range +/- pi
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rc, dc] = iauAticqn(ri, di, astrom, n, b)

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

% Light deflection, giving BCRS coordinate direction.
d = zeros(3,1);
d = iauZp(d);
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
    after = iauLdn(n, b, astrom.eb, before);
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

