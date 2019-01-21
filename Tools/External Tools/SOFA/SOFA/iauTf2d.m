%  - - - - - - - -
%   i a u T f 2 d
%  - - - - - - - -
%
%  Convert hours, minutes, seconds to days.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     s           sign:  '-' = negative, otherwise positive
%     ihour       hours
%     imin        minutes
%     sec         seconds
%
%  Returned:
%     days        interval in days
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
%  This revision:  2012 February 13
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function days = iauTf2d(s, ihour, imin, sec)

constants

% Compute the interval.
if (s == '-')
    temp = -1;
else
    temp = 1;
end

days = temp*( 60*( 60*abs(ihour)+abs(imin) )+abs(sec) )/DAYSEC;

