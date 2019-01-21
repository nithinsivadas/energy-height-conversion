%  - - - - - - - - - -
%   i a u T a i u t 1
%  - - - - - - - - - -
%
%  Time scale transformation:  International Atomic Time, TAI, to
%  Universal Time, UT1.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tai1,tai2      TAI as a 2-part Julian Date
%     dta            UT1-TAI in seconds
%
%  Returned:
%     ut11,ut12      UT1 as a 2-part Julian Date
%
%  Notes:
%
%  1) tai1+tai2 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where tai1 is the Julian
%     Day Number and tai2 is the fraction of a day.  The returned
%     UT11,UT12 follow suit.
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
function [ut11, ut12] = iauTaiut1(tai1, tai2, dta)

constants

% Result, safeguarding precision.
dtad = dta / DAYSEC;

if ( tai1 > tai2 )
    ut11 = tai1;
    ut12 = tai2 + dtad;
else
    ut11 = tai1 + dtad;
    ut12 = tai2;
end

