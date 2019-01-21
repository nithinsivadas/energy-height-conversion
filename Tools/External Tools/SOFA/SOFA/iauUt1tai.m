%  - - - - - - - - - -
%   i a u U t 1 t a i
%  - - - - - - - - - -
%
%  Time scale transformation:  Universal Time, UT1, to International
%  Atomic Time, TAI.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     ut11,ut12      UT1 as a 2-part Julian Date
%     dta            UT1-TAI in seconds
%
%  Returned:
%     tai1,tai2      TAI as a 2-part Julian Date
%
%  Notes:
%  1) ut11+ut12 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where ut11 is the Julian
%     Day Number and ut12 is the fraction of a day.  The returned
%     tai1,tai2 follow suit.
%
%  2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
%     available from IERS tabulations.
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
function [tai1, tai2] = iauUt1tai(ut11, ut12, dta)

constants

% Result, safeguarding precision.
dtad = dta / DAYSEC;

if (ut11 > ut12)
    tai1 = ut11;
    tai2 = ut12 - dtad;
else
    tai1 = ut11 - dtad;
    tai2 = ut12;
end

