%  - - - - - - - -
%   i a u A f 2 a
%  - - - - - - - -
%
%  Convert degrees, arcminutes, arcseconds to radians.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     s           sign:  '-' = negative, otherwise positive
%     ideg        degrees
%     iamin       arcminutes
%     asec        arcseconds
%
%  Returned:
%     rad         angle in radians
%
%  Notes:
%  1)  The result is computed even if any of the range checks fail.
%
%  2)  Negative ideg, iamin and/or asec produce a warning status, but
%      the absolute value is used in the conversion.
%
%  3)  If there are multiple errors, the status value reflects only the
%      first, the smallest taking precedence.
%
%  This revision:  2012 February 13
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rad = iauAf2a(s, ideg, iamin, asec)

constants

% Compute the interval.
if (s == '-')
    temp = -1;
else
    temp =1;
end

rad = temp*( 60*( 60*abs(ideg) + abs(iamin) ) + abs(asec) )*DAS2R;

