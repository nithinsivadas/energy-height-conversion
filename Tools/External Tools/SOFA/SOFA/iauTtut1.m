%  - - - - - - - - -
%   i a u T t u t 1
%  - - - - - - - - -
%
%  Time scale transformation:  Terrestrial Time, TT, to Universal Time,
%  UT1.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tt1,tt2        TT as a 2-part Julian Date
%     dt             TT-UT1 in seconds
%
%  Returned:
%     ut11,ut12      UT1 as a 2-part Julian Date
%
%  Notes:
%  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
%     the two arguments, for example where tt1 is the Julian Day Number
%     and tt2 is the fraction of a day.  The returned ut11,ut12 follow
%     suit.
%
%  2) The argument dt is classical Delta T.
%
%  Reference:
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992)
%
%  This revision:  2011 May 14
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ut11, ut12] = iauTtut1(tt1, tt2, dt)

constants

% Result, safeguarding precision.
dtd = dt / DAYSEC;

if (tt1 > tt2)
    ut11 = tt1;
    ut12 = tt2 - dtd;
else
    ut11 = tt1 - dtd;
    ut12 = tt2;
end

