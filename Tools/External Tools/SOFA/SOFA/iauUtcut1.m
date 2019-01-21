%  - - - - - - - - - -
%   i a u U t c u t 1
%  - - - - - - - - - -
%
%  Time scale transformation:  Coordinated Universal Time, UTC, to
%  Universal Time, UT1.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     utc1,utc2     UTC as a 2-part quasi Julian Date (Notes 1-4)
%     dut1          Delta UT1 = UT1-UTC in seconds (Note 5)
%
%  Returned:
%     ut11,ut12     UT1 as a 2-part Julian Date (Note 6)
%
%  Notes:
%  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
%     convenient way between the two arguments, for example where utc1
%     is the Julian Day Number and utc2 is the fraction of a day.
%
%  2) JD cannot unambiguously represent UTC during a leap second unless
%     special measures are taken.  The convention in the present
%     function is that the JD day represents UTC days whether the
%     length is 86399, 86400 or 86401 SI seconds.
%
%  3) The warning status "dubious year" flags UTCs that predate the
%     introduction of the time scale and that are too far in the future
%     to be trusted.  See iauDat  for further details.
%
%  4) The function iauDtf2d  converts from calendar date and time of
%     day into 2-part Julian Date, and in the case of UTC implements
%     the leap-second-ambiguity convention described above.
%
%  5) Delta UT1 can be obtained from tabulations provided by the
%     International Earth Rotation and Reference Systems Service.  It
%     It is the caller's responsibility to supply a DUT argument
%     containing the UT1-UTC value that matches the given UTC.
%
%  6) The returned ut11,ut12 are such that their sum is the UT1 Julian
%     Date.
%
%  7) The warning status "dubious year" flags UTCs that predate the
%     introduction of the time scale and that are too far in the future
%     to be trusted.  See iauDat for further details.
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992)
%
%  Called:
%     iauJd2cal    JD to Gregorian calendar
%     iauDat       delta(AT) = TAI-UTC
%     iauUtctai    UTC to TAI
%     iauTaiut1    TAI to UT1
%
%  This revision:  2010 May 16
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ut11, ut12] = iauUtcut1(utc1, utc2, dut1)

% Look up TAI-UTC.
[iy, im, id, w] = iauJd2cal(utc1, utc2);
dat = iauDat (iy, im, id, 0);

% Form UT1-TAI.
dta = dut1 - dat;

% UTC to TAI to UT1.
[tai1, tai2] = iauUtctai(utc1, utc2);

[ut11, ut12] = iauTaiut1(tai1, tai2, dta);

