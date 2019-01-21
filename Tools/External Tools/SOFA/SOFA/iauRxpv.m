%  - - - - - - - -
%   i a u R x p v
%  - - - - - - - -
%
%  Multiply a pv-vector by an r-matrix.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     r(3,3)            r-matrix
%     pv(2,3)           pv-vector
%
%  Returned:
%     rpv(2,3)          r * pv
%
%  Note:
%     It is permissible for pv and rpv to be the same array.
%
%  Called:
%     iauRxp       product of r-matrix and p-vector
%
%  This revision:  2013 June 18
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rpv = iauRxpv(r, pv)

rpv(1,:) = iauRxp(r, pv(1,:));
rpv(2,:) = iauRxp(r, pv(2,:));

