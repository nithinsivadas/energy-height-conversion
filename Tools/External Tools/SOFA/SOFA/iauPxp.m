%  - - - - - - -
%   i a u P x p
%  - - - - - - -
%
%  p-vector outer (=vector=cross) product.
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
%     axb(3)            a x b
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
function axb = iauPxp(a, b)

xa = a(1);
ya = a(2);
za = a(3);
xb = b(1);
yb = b(2);
zb = b(3);
axb(1) = ya*zb - za*yb;
axb(2) = za*xb - xa*zb;
axb(3) = xa*yb - ya*xb;

