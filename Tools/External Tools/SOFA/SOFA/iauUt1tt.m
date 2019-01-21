%  - - - - - - - - -
%   i a u U t 1 t t
%  - - - - - - - - -
%
%  Time scale transformation:  Universal Time, UT1, to Terrestrial
%  Time, TT.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     ut11,ut12      UT1 as a 2-part Julian Date
%     dt             TT-UT1 in seconds
%
%  Returned:
%     tt1,tt2        TT as a 2-part Julian Date
%
%  Notes:
%  1) ut11+ut12 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where ut11 is the Julian
%     Day Number and ut12 is the fraction of a day.  The returned
%     tt1,tt2 follow suit.
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
function [tt1, tt2] = iauUt1tt(ut11, ut12, dt)

constants

% Result, safeguarding precision.
dtd = dt / DAYSEC;

if (ut11 > ut12)
    tt1 = ut11;
    tt2 = ut12 + dtd;
else
    tt1 = ut11 + dtd;
    tt2 = ut12;
end

