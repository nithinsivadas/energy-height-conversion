%  - - - - - - - - - -
%   i a u U t 1 u t c
%  - - - - - - - - - -
%
%  Time scale transformation:  Universal Time, UT1, to Coordinated
%  Universal Time, UTC.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     ut11,ut12     UT1 as a 2-part Julian Date (Note 1)
%     dut1          Delta UT1: UT1-UTC in seconds (Note 2)
%
%  Returned:
%     utc1,utc2     UTC as a 2-part quasi Julian Date (Notes 3,4)
%
%  Notes:
%  1) ut11+ut12 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where ut11 is the Julian
%     Day Number and ut12 is the fraction of a day.  The returned utc1
%     and utc2 form an analogous pair, except that a special convention
%     is used, to deal with the problem of leap seconds - see Note 3.
%
%  2) Delta UT1 can be obtained from tabulations provided by the
%     International Earth Rotation and Reference Systems Service.  The
%     value changes abruptly by 1s at a leap second;  however, close to
%     a leap second the algorithm used here is tolerant of the "wrong"
%     choice of value being made.
%
%  3) JD cannot unambiguously represent UTC during a leap second unless
%     special measures are taken.  The convention in the present
%     function is that the returned quasi JD day UTC1+UTC2 represents
%     UTC days whether the length is 86399, 86400 or 86401 SI seconds.
%
%  4) The function iauD2dtf can be used to transform the UTC quasi-JD
%     into calendar date and clock time, including UTC leap second
%     handling.
%
%  5) The warning status "dubious year" flags UTCs that predate the
%     introduction of the time scale and that are too far in the future
%     to be trusted.  See iauDat for further details.
%
%  Called:
%     iauJd2cal    JD to Gregorian calendar
%     iauDat       delta(AT) = TAI-UTC
%     iauCal2jd    Gregorian calendar to JD
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992)
%
%  This revision:  2011 May 14
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [utc1, utc2] = iauUt1utc(ut11, ut12, dut1)

constants

% UT1-UTC in seconds.
duts = dut1;

% Put the two parts of the UT1 into big-first order.
big1 = (ut11 >= ut12);

if (big1)
    u1 = ut11;
    u2 = ut12;
else
    u1 = ut12;
    u2 = ut11;
end

% See if the UT1 can possibly be in a leap-second day.
d1 = u1;
dats1 = 0;
for i = -1:3
    d2 = u2 + i;
    [iy, im, id, fd] = iauJd2cal(d1, d2);
    dats2 = iauDat(iy, im, id, 0);
    if (i == -1)
        dats1 = dats2;
    end
    ddats = dats2 - dats1;
    if (abs(ddats) >= 0.5)
        
        % Yes, leap second nearby: ensure UT1-UTC is "before" value.
        if (ddats * duts >= 0)
            duts = duts - ddats;
        end
        
        % UT1 for the start of the UTC day that ends in a leap.
        [d1, d2] = iauCal2jd(iy, im, id);
        us1 = d1;
        us2 = d2 - 1.0 + duts/DAYSEC;
        
        % Is the UT1 after this point?
        du = u1 - us1;
        du = du + u2 - us2;
        if (du > 0)
            
            % Yes:  fraction of the current UTC day that has elapsed.
            fd = du * DAYSEC / ( DAYSEC + ddats );
            
            % Ramp UT1-UTC to bring about SOFA's JD(UTC) convention.
            if (fd <= 1)
                temp = fd;
            else
                temp =1;
            end
            
            duts = duts + ddats * temp;
        end
        
        % Done.
        break
    end
    dats1 = dats2;
end

% Subtract the (possibly adjusted) UT1-UTC from UT1 to give UTC.
u2 = u2 - duts / DAYSEC;

% Result, safeguarding precision.
if (big1)
    utc1 = u1;
    utc2 = u2;
else
    utc1 = u2;
    utc2 = u1;
end

