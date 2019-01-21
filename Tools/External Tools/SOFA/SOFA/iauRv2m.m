%  - - - - - - - -
%   i a u R v 2 m
%  - - - - - - - -
%
%  Form the r-matrix corresponding to a given r-vector.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     w(3)              rotation vector (Note 1)
%
%  Returned:
%     r(3,3)            rotation matrix
%
%  Notes:
%  1) A rotation matrix describes a rotation through some angle about
%     some arbitrary axis called the Euler axis.  The "rotation vector"
%     supplied to This function has the same direction as the Euler
%     axis, and its magnitude is the angle in radians.
%
%  2) If w is null, the unit matrix is returned.
%
%  3) The reference frame rotates clockwise as seen looking along the
%     rotation vector from the origin.
%
%  This revision:  2008 May 11
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = iauRv2m(w)

% Euler angle (magnitude of rotation vector) and functions.
x = w(1);
y = w(2);
z = w(3);
phi = sqrt(x*x + y*y + z*z);
s = sin(phi);
c = cos(phi);
f = 1 - c;

% Euler axis (direction of rotation vector), perhaps null.
if (phi ~= 0)
    x = x/phi;
    y = y/phi;
    z = z/phi;
end

% Form the rotation matrix.
r(1,1) = x*x*f + c;
r(1,2) = x*y*f + z*s;
r(1,3) = x*z*f - y*s;
r(2,1) = y*x*f - z*s;
r(2,2) = y*y*f + c;
r(2,3) = y*z*f + x*s;
r(3,1) = z*x*f + y*s;
r(3,2) = z*y*f - x*s;
r(3,3) = z*z*f + c;

