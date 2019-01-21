%  - - - - - - -
%   i a u S 2 c
%  - - - - - - -
%
%  Convert spherical coordinates to Cartesian.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     theta           longitude angle (radians)
%     phi             latitude angle (radians)
%
%  Returned:
%     c(3)            direction cosines
%
%  This revision:  2008 October 28
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = iauS2c(theta, phi)

cp = cos(phi);
c(1) = cos(theta) * cp;
c(2) = sin(theta) * cp;
c(3) = sin(phi);

