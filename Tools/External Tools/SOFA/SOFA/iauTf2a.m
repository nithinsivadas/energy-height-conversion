%  - - - - - - - -
%   i a u T f 2 a
%  - - - - - - - -
%
%  Convert hours, minutes, seconds to radians.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     s         char    sign:  '-' = negative, otherwise positive
%     ihour     int     hours
%     imin      int     minutes
%     sec       double  seconds
%
%  Returned:
%     rad       double  angle in radians
%
%  Notes:
%  1)  The result is computed even if any of the range checks fail.
%
%  2)  Negative ihour, imin and/or sec produce a warning status, but
%      the absolute value is used in the conversion.
%
%  3)  If there are multiple errors, the status value reflects only the
%      first, the smallest taking precedence.
%
%  This revision:  2013 June 18
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rad = iauTf2a(s, ihour, imin, sec)

constants

% Compute the interval.
if ( s == '-')
    temp = -1;
else
    temp = 1;
end

rad  = temp*( 60*( 60*( abs(ihour) )+( abs(imin) ) ) + abs(sec) )*DS2R;

