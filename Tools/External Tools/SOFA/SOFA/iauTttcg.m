%  - - - - - - - - -
%   i a u T t t c g
%  - - - - - - - - -
%
%  Time scale transformation:  Terrestrial Time, TT, to Geocentric
%  Coordinate Time, TCG.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tt1,tt2        TT as a 2-part Julian Date
%
%  Returned:
%     tcg1,tcg2      TCG as a 2-part Julian Date
%
%  Note:
%     tt1+tt2 is Julian Date, apportioned in any convenient way between
%     the two arguments, for example where tt1 is the Julian Day Number
%     and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
%     suit.
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%     IAU 2000 Resolution B1.9
%
%  This revision:  2010 May 13
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tcg1, tcg2] = iauTttcg(tt1, tt2)

constants

% 1977 Jan 1 00:00:32.184 TT, as MJD
t77t = DJM77 + TTMTAI/DAYSEC;

% TT to TCG rate
elgg = ELG/(1-ELG);

% Result, safeguarding precision.
if (tt1 > tt2)
    tcg1 = tt1;
    tcg2 = tt2 + ( ( tt1 - DJM0 ) + ( tt2 - t77t ) ) * elgg;
else
    tcg1 = tt1 + ( ( tt2 - DJM0 ) + ( tt1 - t77t ) ) * elgg;
    tcg2 = tt2;
end

