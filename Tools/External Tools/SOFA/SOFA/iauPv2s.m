%  - - - - - - - -
%   i a u P v 2 s
%  - - - - - - - -
%
%  Convert position/velocity from Cartesian to spherical coordinates.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     pv(2,3)          pv-vector
%
%  Returned:
%     theta            longitude angle (radians)
%     phi              latitude angle (radians)
%     r                radial distance
%     td               rate of change of theta
%     pd               rate of change of phi
%     rd               rate of change of r
%
%  Notes:
%  1) If the position part of pv is null, theta, phi, td and pd
%     are indeterminate.  This is handled by extrapolating the
%     position through unit time by using the velocity part of
%     pv.  This moves the origin without changing the direction
%     of the velocity component.  If the position and velocity
%     components of pv are both null, zeroes are returned for all
%     six results.
%
%  2) If the position is a pole, theta, td and pd are indeterminate.
%     In such cases zeroes are returned for all three.
%
%  This revision:  2008 October 28
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, phi, r, td, pd, rd] = iauPv2s(pv)

% Components of position/velocity vector.
x  = pv(1,1);
y  = pv(1,2);
z  = pv(1,3);
xd = pv(2,1);
yd = pv(2,2);
zd = pv(2,3);

% Component of r in XY plane squared.
rxy2 = x*x + y*y;

% Modulus squared.
r2 = rxy2 + z*z;

% Modulus.
rtrue = sqrt(r2);

% If null vector, move the origin along the direction of movement.
rw = rtrue;
if (rtrue == 0)
    x = xd;
    y = yd;
    z = zd;
    rxy2 = x*x + y*y;
    r2 = rxy2 + z*z;
    rw = sqrt(r2);
end

% Position and velocity in spherical coordinates.
rxy = sqrt(rxy2);
xyp = x*xd + y*yd;
if (rxy2 ~= 0)
    theta = atan2(y, x);
    phi = atan2(z, rxy);
    td = (x*yd - y*xd) / rxy2;
    pd = (zd*rxy2 - z*xyp) / (r2*rxy);
else
    theta = 0;
    if (z ~= 0)
        phi = atan2(z, rxy);
    else
        phi = 0;
    end
    td = 0;
    pd = 0;
end
r = rtrue;
if (rw ~= 0)
    rd = (xyp + z*zd) / rw;
else
    rd = 0;
end

