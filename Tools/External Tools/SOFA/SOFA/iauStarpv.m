%  - - - - - - - - - -
%   i a u S t a r p v
%  - - - - - - - - - -
%
%  Convert star catalog coordinates to position+velocity vector.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given (Note 1):
%     ra             right ascension (radians)
%     dec            declination (radians)
%     pmr            RA proper motion (radians/year)
%     pmd            Dec proper motion (radians/year)
%     px             parallax (arcseconds)
%     rv             radial velocity (km/s, positive = receding)
%
%  Returned (Note 2):
%     pv(2,3)        pv-vector (AU, AU/day)
%
%  Notes:
%  1) The star data accepted by this function are "observables" for an
%     imaginary observer at the solar-system barycenter.  Proper motion
%     and radial velocity are, strictly, in terms of barycentric
%     coordinate time, TCB.  For most practical applications, it is
%     permissible to neglect the distinction between TCB and ordinary
%     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
%     limited by the intrinsic accuracy of the proper-motion and
%     radial-velocity data;  moreover, the pv-vector is likely to be
%     merely an intermediate result, so that a change of time unit
%     would cancel out overall.
%
%     In accordance with normal star-catalog conventions, the object's
%     right ascension and declination are freed from the effects of
%     secular aberration.  The frame, which is aligned to the catalog
%     equator and equinox, is Lorentzian and centered on the SSB.
%
%  2) The resulting position and velocity pv-vector is with respect to
%     the same frame and, like the catalog coordinates, is freed from
%     the effects of secular aberration.  Should the "coordinate
%     direction", where the object was located at the catalog epoch, be
%     required, it may be obtained by calculating the magnitude of the
%     position vector pv[0][0-2] dividing by the speed of light in
%     AU/day to give the light-time, and then multiplying the space
%     velocity pv[1][0-2] by this light-time and adding the result to
%     pv[0][0-2].
%
%     Summarizing, the pv-vector returned is for most stars almost
%     identical to the result of applying the standard geometrical
%     "space motion" transformation.  The differences, which are the
%     subject of the Stumpff paper referenced below, are:
%
%     (i) In stars with significant radial velocity and proper motion,
%     the constantly changing light-time distorts the apparent proper
%     motion.  Note that this is a classical, not a relativistic,
%     effect.
%
%     (ii) The transformation complies with special relativity.
%
%  3) Care is needed with units.  The star coordinates are in radians
%     and the proper motions in radians per Julian year, but the
%     parallax is in arcseconds; the radial velocity is in km/s, but
%     the pv-vector result is in AU and AU/day.
%
%  4) The RA proper motion is in terms of coordinate angle, not true
%     angle.  If the catalog uses arcseconds for both RA and Dec proper
%     motions, the RA proper motion will need to be divided by cos(Dec)
%     before use.
%
%  5) Straight-line motion at constant speed, in the inertial frame,
%     is assumed.
%
%  6) An extremely small (or zero or negative) parallax is interpreted
%     to mean that the object is on the "celestial sphere", the radius
%     of which is an arbitrary (large) value (see the constant PXMIN).
%     When the distance is overridden in this way, the status,
%     initially zero, has 1 added to it.
%
%  7) If the space velocity is a significant fraction of c (see the
%     constant VMAX), it is arbitrarily set to zero.  When this action
%     occurs, 2 is added to the status.
%
%  8) The relativistic adjustment involves an iterative calculation.
%     If the process fails to converge within a set number (IMAX) of
%     iterations, 4 is added to the status.
%
%  9) The inverse transformation is performed by the function
%     iauPvstar.
%
%  Called:
%     iauS2pv      spherical coordinates to pv-vector
%     iauPm        modulus of p-vector
%     iauZp        zero p-vector
%     iauPn        decompose p-vector into modulus and direction
%     iauPdp       scalar product of two p-vectors
%     iauSxp       multiply p-vector by scalar
%     iauPmp       p-vector minus p-vector
%     iauPpp       p-vector plus p-vector
%
%  Reference:
%     Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.
%
%  This revision:  2009 July 6
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pv = iauStarpv(ra, dec, pmr, pmd, px, rv)

constants

% Smallest allowed parallax
PXMIN = 1e-7;

% Largest allowed speed (fraction of c)
VMAX = 0.5;

% Maximum number of iterations for relativistic solution
IMAX = 100;

d = 0;
del = 0;
odd = 0;
oddel = 0;
od = 0;
odel = 0;

% Distance (AU).
if (px >= PXMIN)
    w = px;    
else
    w = PXMIN;    
end
r = DR2AS / w;

% Radial velocity (AU/day).
rd = DAYSEC * rv * 1e3 / DAU;

% Proper motion (radian/day).
rad = pmr / DJY;
decd = pmd / DJY;

% To pv-vector (AU,AU/day).
pv = iauS2pv(ra, dec, r, rad, decd, rd);

% If excessive velocity, arbitrarily set it to zero.
v = iauPm(pv(2,:));
if (v / DC > VMAX)
    pv(2,:) = iauZp(pv(2,:));    
end

% Isolate the radial component of the velocity (AU/day).
[w, x] = iauPn(pv(1,:));
vsr = iauPdp(x, pv(2,:));
usr = iauSxp(vsr, x);

% Isolate the transverse component of the velocity (AU/day).
ust = iauPmp(pv(2,:), usr);
vst = iauPm(ust);

% Special-relativity dimensionless parameters.
betsr = vsr / DC;
betst = vst / DC;

% Determine the inertial-to-observed relativistic correction terms.
bett = betst;
betr = betsr;
for i = 0:IMAX-1
    d = 1 + betr;
    del = sqrt(1 - betr*betr - bett*bett) - 1;
    betr = d * betsr + del;
    bett = d * betst;
    if (i > 0)
        dd = abs(d - od);
        ddel = abs(del - odel);
        if ((i > 1) && (dd >= odd) && (ddel >= oddel))
            break
        end
        odd = dd;
        oddel = ddel;
    end
    od = d;
    odel = del;
end

% Replace observed radial velocity with inertial value.
if (betsr ~= 0)
    w = d + del / betsr;
else
    w = 1;
end

ur = iauSxp(w, usr);

% Replace observed tangential velocity with inertial value.
ut = iauSxp(d, ust);

% Combine the two to obtain the inertial space velocity.
pv(2,:) = iauPpp(ur, ut);

