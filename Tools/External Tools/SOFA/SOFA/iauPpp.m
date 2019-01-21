%  - - - - - - -
%   i a u P p p
%  - - - - - - -
%
%  P-vector addition.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     a(3)              first p-vector
%     b(3)              second p-vector
%
%  Returned:
%     apb(3)            a + b
%
%  Note:
%     It is permissible to re-use the same array for any of the
%     arguments.
%
%  This revision:  2008 November 18
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function apb = iauPpp(a, b)

apb(1) = a(1) + b(1);
apb(2) = a(2) + b(2);
apb(3) = a(3) + b(3);

