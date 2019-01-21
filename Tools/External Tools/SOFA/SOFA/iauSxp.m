%  - - - - - - -
%   i a u S x p
%  - - - - - - -
%
%  Multiply a p-vector by a scalar.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     s              scalar
%     p(3)           p-vector
%
%  Returned:
%     sp(3)          s * p
%
%  Note:
%     It is permissible for p and sp to be the same array.
%
%  This revision:  2008 October 28
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sp = iauSxp(s, p)

sp(1) = s * p(1);
sp(2) = s * p(2);
sp(3) = s * p(3);

