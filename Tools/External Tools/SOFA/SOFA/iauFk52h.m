%  - - - - - - - - -
%   i a u F k 5 2 h
%  - - - - - - - - -
%
%  Transform FK5 (J2000.0) star data into the Hipparcos system.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given (all FK5, equinox J2000.0, epoch J2000.0):
%     r5         RA (radians)
%     d5         Dec (radians)
%     dr5        proper motion in RA (dRA/dt, rad/Jyear)
%     dd5        proper motion in Dec (dDec/dt, rad/Jyear)
%     px5        parallax (arcsec)
%     rv5        radial velocity (km/s, positive = receding)
%
%  Returned (all Hipparcos, epoch J2000.0):
%     rh         RA (radians)
%     dh         Dec (radians)
%     drh        proper motion in RA (dRA/dt, rad/Jyear)
%     ddh        proper motion in Dec (dDec/dt, rad/Jyear)
%     pxh        parallax (arcsec)
%     rvh        radial velocity (km/s, positive = receding)
%
%  Notes:
%  1) This function transforms FK5 star positions and proper motions
%     into the system of the Hipparcos catalog.
%
%  2) The proper motions in RA are dRA/dt rather than
%     cos(Dec)*dRA/dt, and are per year rather than per century.
%
%  3) The FK5 to Hipparcos transformation is modeled as a pure
%     rotation and spin;  zonal errors in the FK5 catalog are not
%     taken into account.
%
%  4) See also iauH2fk5, iauFk5hz, iauHfk5z.
%
%  Called:
%     iauStarpv    star catalog data to space motion pv-vector
%     iauFk5hip    FK5 to Hipparcos rotation and spin
%     iauRxp       product of r-matrix and p-vector
%     iauPxp       vector product of two p-vectors
%     iauPpp       p-vector plus p-vector
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
function [rh, dh, drh, ddh, pxh, rvh] = iauFk52h(r5, d5, dr5, dd5, px5, rv5)

% FK5 barycentric position/velocity pv-vector (normalized).
pv5 = iauStarpv(r5, d5, dr5, dd5, px5, rv5);

% FK5 to Hipparcos orientation matrix and spin vector.
r5h = zeros(3);
s5h = zeros(3,1);
[r5h, s5h] = iauFk5hip(r5h, s5h);

% Make spin units per day instead of per year.
for i = 1:3
    s5h(i) = s5h(i)/365.25;
end

% Orient the FK5 position into the Hipparcos system.
pvh(1,:) = iauRxp(r5h, pv5(1,:));

% Apply spin to the position giving an extra space motion component.
wxp = iauPxp(pv5(1,:), s5h);

% Add this component to the FK5 space motion.
vv = iauPpp(wxp, pv5(2,:));

% Orient the FK5 space motion into the Hipparcos system.
pvh(2,:) = iauRxp(r5h, vv);

% Hipparcos pv-vector to spherical.
[rh, dh, drh, ddh, pxh, rvh] = iauPvstar(pvh);

