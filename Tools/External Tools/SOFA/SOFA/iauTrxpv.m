%  - - - - - - - - -
%   i a u T r x p v
%  - - - - - - - - -
%
%  Multiply a pv-vector by the transpose of an r-matrix.
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
%     trpv(2,3)         r * pv
%
%  Note:
%     It is permissible for pv and trpv to be the same array.
%
%  Called:
%     iauTr        transpose r-matrix
%     iauRxpv      product of r-matrix and pv-vector
%
%  This revision:  2013 June 18
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trpv = iauTrxpv(r, pv)

% Transpose of matrix r.
tr = iauTr(r);

% Matrix tr * vector pv -> vector trpv.
trpv = iauRxpv(tr, pv);

