%  - - - - - - -
%   i a u C 2 s
%  - - - - - - -
%
%  P-vector to spherical coordinates.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     p(3)          p-vector
%
%  Returned:
%     theta         longitude angle (radians)
%     phi           latitude angle (radians)
%
%  Notes:
%
%  1) The vector p can have any magnitude; only its direction is used.
%
%  2) If p is null, zero theta and phi are returned.
%
%  3) At either pole, zero theta is returned.
%
%  This revision:  2008 May 11
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, phi] = iauC2s(p)

x  = p(1);
y  = p(2);
z  = p(3);
d2 = x*x + y*y;

if(d2 == 0)
    theta = 0;
else
    theta = atan2(y, x);
end
if(z == 0)
    phi = 0;
else
    phi = atan2(z, sqrt(d2));
end

