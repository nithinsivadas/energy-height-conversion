%  - - - - - - -
%   i a u P m p
%  - - - - - - -
%
%  P-vector subtraction.
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
%     amb(3)            a - b
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
function amb = iauPmp(a, b)

amb(1) = a(1) - b(1);
amb(2) = a(2) - b(2);
amb(3) = a(3) - b(3);

