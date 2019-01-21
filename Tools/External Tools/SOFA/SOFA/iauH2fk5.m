%  - - - - - - - - -
%   i a u H 2 f k 5
%  - - - - - - - - -
%
%  Transform Hipparcos star data into the FK5 (J2000.0) system.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given (all Hipparcos, epoch J2000.0):
%     rh         RA (radians)
%     dh         Dec (radians)
%     drh        proper motion in RA (dRA/dt, rad/Jyear)
%     ddh        proper motion in Dec (dDec/dt, rad/Jyear)
%     pxh        parallax (arcsec)
%     rvh        radial velocity (km/s, positive = receding)
%
%  Returned (all FK5, equinox J2000.0, epoch J2000.0):
%     r5         RA (radians)
%     d5         Dec (radians)
%     dr5        proper motion in RA (dRA/dt, rad/Jyear)
%     dd5        proper motion in Dec (dDec/dt, rad/Jyear)
%     px5        parallax (arcsec)
%     rv5        radial velocity (km/s, positive = receding)
%
%  Notes:
%  1) This function transforms Hipparcos star positions and proper
%     motions into FK5 J2000.0.
%
%  2) The proper motions in RA are dRA/dt rather than
%     cos(Dec)*dRA/dt, and are per year rather than per century.
%
%  3) The FK5 to Hipparcos transformation is modeled as a pure
%     rotation and spin;  zonal errors in the FK5 catalog are not
%     taken into account.
%
%  4) See also iauFk52h, iauFk5hz, iauHfk5z.
%
%  Called:
%     iauStarpv    star catalog data to space motion pv-vector
%     iauFk5hip    FK5 to Hipparcos rotation and spin
%     iauRv2m      r-vector to r-matrix
%     iauRxp       product of r-matrix and p-vector
%     iauTrxp      product of transpose of r-matrix and p-vector
%     iauPxp       vector product of two p-vectors
%     iauPmp       p-vector minus p-vector
%     iauPvstar    space motion pv-vector to star catalog data
%
%  Reference:
%     F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).
%
%  This revision:  2009 December 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r5, d5, dr5, dd5, px5, rv5] = iauH2fk5(rh, dh, drh, ddh, pxh, rvh)

% Hipparcos barycentric position/velocity pv-vector (normalized).
pvh = iauStarpv(rh, dh, drh, ddh, pxh, rvh);

% FK5 to Hipparcos orientation matrix and spin vector.
r5h = zeros(3);
s5h = zeros(3,1);
[r5h, s5h] = iauFk5hip(r5h, s5h);

% Make spin units per day instead of per year.
for i = 1:3
    s5h(i) = s5h(i)/365.25;
end

% Orient the spin into the Hipparcos system.
sh = iauRxp(r5h, s5h);

% De-orient the Hipparcos position into the FK5 system.
pv5(1,:) = iauTrxp(r5h, pvh(1,:));

% Apply spin to the position giving an extra space motion component.
wxp = iauPxp(pvh(1,:), sh);

% Subtract this component from the Hipparcos space motion.
vv = iauPmp(pvh(1,:), wxp);

% De-orient the Hipparcos space motion into the FK5 system.
pv5(2,:) = iauTrxp(r5h, vv);

% FK5 pv-vector to spherical.
[r5, d5, dr5, dd5, px5, rv5] = iauPvstar(pv5);

