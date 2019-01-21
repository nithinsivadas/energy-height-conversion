%  - - - - - -
%   i a u P n
%  - - - - - -
%
%  Convert a p-vector into modulus and unit vector.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     p(3)              p-vector
%
%  Returned:
%     r                 modulus
%     u(3)              unit vector
%
%  Notes:
%  1) If p is null, the result is null.  Otherwise the result is a unit
%     vector.
%
%  2) It is permissible to re-use the same array for any of the
%     arguments.
%
%  Called:
%     iauPm        modulus of p-vector
%     iauZp        zero p-vector
%     iauSxp       multiply p-vector by scalar
%
%  This revision:  2008 November 18
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, u] = iauPn(p)

% Obtain the modulus and test for zero.
w = iauPm(p);

if (w == 0)
    % Null vector.
    u = zeros(3,1);
    u = iauZp(u);
else
    % Unit vector.
    u = iauSxp(1/w, p);
end

% Return the modulus.
r = w;

