%  - - - - - - -
%   i a u P d p
%  - - - - - - -
%
%  p-vector inner (=scalar=dot) product.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     a(3)           first p-vector
%     b(3)           second p-vector
%
%  Returned (function value):
%                    a . b
%
%  This revision:  2008 May 22
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = iauPdp(a, b)

w  = a(1) * b(1) + a(2) * b(2) + a(3) * b(3);

