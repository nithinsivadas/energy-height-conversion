%  - - - - - - - - -
%   i a u A t o i q
%  - - - - - - - - -
%
%  Quick observed place to CIRS, given the star-independent astrometry
%  parameters.
%
%  Use of this function is appropriate when efficiency is important and
%  where many star positions are all to be transformed for one date.
%  The star-independent astrometry parameters can be obtained by
%  calling iauApio[13] or iauApco[13].
%
%  Status:  support function.
%
%  Given:
%     type        type of coordinates: "R", "H" or "A" (Note 1)
%     ob1         observed Az, HA or RA (radians; Az is N=0,E=90)
%     ob2         observed ZD or Dec (radians)
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
%     ri          CIRS right ascension (CIO-based, radians)
%     di          CIRS declination (radians)
%
%  Notes:
%  1) "Observed" Az,El means the position that would be seen by a
%     perfect geodetically aligned theodolite.  This is related to
%     the observed HA,Dec via the standard rotation, using the geodetic
%     latitude (corrected for polar motion), while the observed HA and
%     RA are related simply through the Earth rotation angle and the
%     site longitude.  "Observed" RA,Dec or HA,Dec thus means the
%     position that would be seen by a perfect equatorial with its
%     polar axis aligned to the Earth's axis of rotation.  By removing
%     from the observed place the effects of atmospheric refraction and
%     diurnal aberration, the CIRS RA,Dec is obtained.
%
%  2) Only the first character of the type argument is significant.
%     "R" or "r" indicates that ob1 and ob2 are the observed right
%     ascension and declination;  "H" or "h" indicates that they are
%     hour angle (west +ve) and declination;  anything else ("A" or
%     "a" is recommended) indicates that ob1 and ob2 are azimuth (north
%     zero, east 90 deg) and zenith distance.  (Zenith distance is used
%     rather than altitude in order to reflect the fact that no
%     allowance is made for depression of the horizon.)
%
%  3) The accuracy of the result is limited by the corrections for
%     refraction, which use a simple A*tan(z) + B*tan^3(z) model.
%     Providing the meteorological parameters are known accurately and
%     there are no gross local effects, the predicted observed
%     coordinates should be within 0.05 arcsec (optical) or 1 arcsec
%     (radio) for a zenith distance of less than 70 degrees, better
%     than 30 arcsec (optical or radio) at 85 degrees and better than
%     20 arcmin (optical) or 30 arcmin (radio) at the horizon.
%
%     Without refraction, the complementary functions iauAtioq and
%     iauAtoiq are self-consistent to better than 1 microarcsecond all
%     over the celestial sphere.  With refraction included, consistency
%     falls off at high zenith distances, but is still better than
%     0.05 arcsec at 85 degrees.
%
%  4) It is advisable to take great care with units, as even unlikely
%     values of the input parameters are accepted and processed in
%     accordance with the models used.
%
%  Called:
%     iauS2c       spherical coordinates to unit vector
%     iauC2s       p-vector to spherical
%     iauAnp       normalize angle into range 0 to 2pi
%
%  This revision:   2013 October 9
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ri, di] = iauAtoiq(type, ob1, ob2, astrom)

% Coordinate type.
c = type;

% Coordinates.
c1 = ob1;
c2 = ob2;

% Sin, cos of latitude.
sphi = astrom.sphi;
cphi = astrom.cphi;

% Standardize coordinate type.
if (c == 'r' || c == 'R')
      c = 'R';
elseif (c == 'h' || c == 'H')
    c = 'H';
else
    c = 'A';
end

% If Az,ZD, convert to Cartesian (S=0,E=90).
if (c == 'A')
    ce = sin(c2);
    xaeo = - cos(c1) * ce;
    yaeo = sin(c1) * ce;
    zaeo = cos(c2);
else
    % If RA,Dec, convert to HA,Dec.
    if (c == 'R')
        c1 = astrom.eral - c1;
    end
    % To Cartesian -HA,Dec.
    v = iauS2c(-c1, c2);
    xmhdo = v(1);
    ymhdo = v(2);
    zmhdo = v(3);    
    % To Cartesian Az,El (S=0,E=90).
    xaeo = sphi*xmhdo - cphi*zmhdo;
    yaeo = ymhdo;
    zaeo = cphi*xmhdo + sphi*zmhdo;
end

% Azimuth (S=0,E=90).
if (xaeo ~= 0 || yaeo ~= 0)
    az = atan2(yaeo, xaeo);
else
    az = 0;
end

% Sine of observed ZD, and observed ZD.
sz = sqrt(xaeo*xaeo + yaeo*yaeo);
zdo = atan2(sz, zaeo);

% Refraction
% ----------

% Fast algorithm using two constant model.
refa = astrom.refa;
refb = astrom.refb;
tz = sz / zaeo;
dref = (refa + refb*tz*tz) * tz;
zdt = zdo + dref;

% To Cartesian Az,ZD.
ce = sin(zdt);
xaet = cos(az) * ce;
yaet = sin(az) * ce;
zaet = cos(zdt);

% Cartesian Az,ZD to Cartesian -HA,Dec.
xmhda = sphi*xaet + cphi*zaet;
ymhda = yaet;
zmhda = - cphi*xaet + sphi*zaet;

% Diurnal aberration.
f = ( 1 + astrom.diurab*ymhda );
xhd = f * xmhda;
yhd = f * ( ymhda - astrom.diurab );
zhd = f * zmhda;

% Polar motion.
xpl = astrom.xpl;
ypl = astrom.ypl;
w = xpl*xhd - ypl*yhd + zhd;
v(1) = xhd - xpl*w;
v(2) = yhd + ypl*w;
v(3) = w - ( xpl*xpl + ypl*ypl ) * zhd;

% To spherical -HA,Dec.
[hma, di] = iauC2s(v);

% Right ascension.
ri = iauAnp(astrom.eral + hma);

