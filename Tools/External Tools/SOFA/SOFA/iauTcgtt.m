%  - - - - - - - - -
%   i a u T c g t t
%  - - - - - - - - -
%
%  Time scale transformation:  Geocentric Coordinate Time, TCG, to
%  Terrestrial Time, TT.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tcg1,tcg2      TCG as a 2-part Julian Date
%
%  Returned:
%     tt1,tt2        TT as a 2-part Julian Date
%
%  Note:
%
%     tcg1+tcg2 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where tcg1 is the Julian
%     Day Number and tcg22 is the fraction of a day.  The returned
%     tt1,tt2 follow suit.
%
%  References:
%
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),.
%     IERS Technical Note No. 32, BKG (2004)
%
%     IAU 2000 Resolution B1.9
%
%  This revision:  2010 May 14
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tt1, tt2] = iauTcgtt(tcg1, tcg2)

constants

% 1977 Jan 1 00:00:32.184 TT, as MJD
t77t = DJM77 + TTMTAI/DAYSEC;

% Result, safeguarding precision.
if (tcg1 > tcg2)
    tt1 = tcg1;
    tt2 = tcg2 - ( ( tcg1 - DJM0 ) + ( tcg2 - t77t ) ) * ELG;
else
    tt1 = tcg1 - ( ( tcg2 - DJM0 ) + ( tcg1 - t77t ) ) * ELG;
    tt2 = tcg2;
end

