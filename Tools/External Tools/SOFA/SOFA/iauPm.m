%  - - - - - -
%   i a u P m
%  - - - - - -
%
%  Modulus of p-vector.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     p(3)           p-vector
%
%  This revision:  2008 May 22
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = iauPm(p)

w  = sqrt( p(1) * p(1) + p(2) * p(2) + p(3) * p(3) );

