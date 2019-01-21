%  - - - - - - - -
%   i a u P p s p
%  - - - - - - - -
%
%  P-vector plus scaled p-vector.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     a(3)           first p-vector
%     s              scalar (multiplier for b)
%     b(3)           second p-vector
%
%  Returned:
%     apsb(3)        a + s*b
%
%  Note:
%     It is permissible for any of a, b and apsb to be the same array.
%
%  Called:
%     iauSxp       multiply p-vector by scalar
%     iauPpp       p-vector plus p-vector
%
%  This revision:  2008 November 18
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function apsb = iauPpsp(a, s, b)

% s*b.
sb = iauSxp(s, b);

% a + s*b.
apsb = iauPpp(a, sb);

