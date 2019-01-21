%  - - - - - - - -
%   i a u T r x p
%  - - - - - - - -
%
%  Multiply a p-vector by the transpose of an r-matrix.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     r(3,3)           r-matrix
%     p(3)             p-vector
%
%  Returned:
%     trp(3)           r * p
%
%  Note:
%     It is permissible for p and trp to be the same array.
%
%  Called:
%     iauTr        transpose r-matrix
%     iauRxp       product of r-matrix and p-vector
%
%  This revision:  2008 October 28
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trp = iauTrxp(r, p)

% Transpose of matrix r.
tr = iauTr(r);

% Matrix tr * vector p -> vector trp.
trp = iauRxp(tr, p);

