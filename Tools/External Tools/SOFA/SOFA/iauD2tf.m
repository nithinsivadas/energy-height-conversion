%  - - - - - - - -
%   i a u D 2 t f
%  - - - - - - - -
%
%  Decompose days to hours, minutes, seconds, fraction.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     ndp        resolution (Note 1)
%     days       interval in days
%
%  Returned:
%     sign       '+' or '-'
%     ihmsf(4)   hours, minutes, seconds, fraction
%
%  1) The argument ndp is interpreted as follows:
%
%     ndp         resolution
%      :      ...0000 00 00
%     -7         1000 00 00
%     -6          100 00 00
%     -5           10 00 00
%     -4            1 00 00
%     -3            0 10 00
%     -2            0 01 00
%     -1            0 00 10
%      0            0 00 01
%      1            0 00 00.1
%      2            0 00 00.01
%      3            0 00 00.001
%      :            0 00 00.000...
%
%  2) The largest positive useful value for ndp is determined by the
%     size of days, the format of double on the target platform, and
%     the risk of overflowing ihmsf[3].  On a typical platform, for
%     days up to 1.0, the available floating-point precision might
%     correspond to ndp=12.  However, the practical limit is typically
%     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
%     only 16 bits.
%
%  3) The absolute value of days may exceed 1.0.  In cases where it
%     does not, it is up to the caller to test for and handle the
%     case where days is very nearly 1.0 and rounds up to 24 hours,
%     by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
%
%  This revision:  2008 May 17
%
%  SOFA release 2012-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sign, ihmsf] = iauD2tf(ndp, days)

constants

% Handle sign.
if ( days >= 0 )
    sign = '+';
else
    sign = '-';
end

% Interval in seconds.
a = DAYSEC * abs(days);

% Pre-round if resolution coarser than 1s (then pretend ndp=1).
if (ndp < 0)
    nrs = 1;
    for n = 1:-ndp
        if (n == 2 || n == 4)
            temp = 6;
        else
            temp = 10;
        end
        nrs = nrs * temp;
    end
    rs = nrs;
    w = a / rs;
    a = rs * dnint(w);
end

% Express the unit of each field in resolution units.
nrs = 1;

for n = 1:ndp
    nrs = nrs * 10;
end

rs = nrs;
rm = rs * 60;
rh = rm * 60;

% Round the interval and express in resolution units.
a = dnint(rs * a);

% Break into fields.
ah = a / rh;
ah = dint(ah);
a = a - ah * rh;
am = a / rm;
am = dint(am);
a = a - am * rm;
as = a / rs;
as = dint(as);
af = a - as * rs;

% Return results.
ihmsf(1) = int64(ah);
ihmsf(2) = int64(am);
ihmsf(3) = int64(as);
ihmsf(4) = int64(af);

