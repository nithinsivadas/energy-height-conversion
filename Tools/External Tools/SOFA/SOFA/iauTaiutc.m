%  - - - - - - - - - -
%   i a u T a i u t c
%  - - - - - - - - - -
%
%  Time scale transformation:  International Atomic Time, TAI, to
%  Coordinated Universal Time, UTC.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tai1,tai2     TAI as a 2-part Julian Date (Note 1)
%
%  Returned:
%     utc1,utc2     UTC as a 2-part quasi Julian Date (Notes 1-3)
%
%  Notes:
%
%  1) tai1+tai2 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where tai1 is the Julian
%     Day Number and tai2 is the fraction of a day.  The returned utc1
%     and utc2 form an analogous pair, except that a special convention
%     is used, to deal with the problem of leap seconds - see the next
%     note.
%
%  2) JD cannot unambiguously represent UTC during a leap second unless
%     special measures are taken.  The convention in the present
%     function is that the JD day represents UTC days whether the
%     length is 86399, 86400 or 86401 SI seconds.
%
%  3) The function iauD2dtf can be used to transform the UTC quasi-JD
%     into calendar date and clock time, including UTC leap second
%     handling.
%
%  4) The warning status "dubious year" flags UTCs that predate the
%     introduction of the time scale and that are too far in the future
%     to be trusted.  See iauDat for further details.
%
%  Called:
%     iauJd2cal    JD to Gregorian calendar
%     iauDat       delta(AT) = TAI-UTC
%     iauCal2jd    Gregorian calendar to JD
%
%  References:
%
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992)
%
%  This revision:  2011 May 14
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [utc1, utc2] = iauTaiutc(tai1, tai2)

constants

% Put the two parts of the TAI into big-first order.
big1 = ( tai1 >= tai2 );

if (big1)
    a1 = tai1;
    a2 = tai2;
else
    a1 = tai2;
    a2 = tai1;
end

% See if the TAI can possibly be in a leap-second day.
d1 = a1;
dats1 = 0;
for i = -1:3
    d2 = a2 + i;
    [iy, im, id, fd] = iauJd2cal(d1, d2);
    dats2 = iauDat(iy, im, id, 0);
    if (i == -1)
        dats1 = dats2;
    end
    ddats = dats2 - dats1;
    datd = dats1 / DAYSEC;
    if (abs(ddats) >= 0.5)
        % Yes.  Get TAI for the start of the UTC day that
        % ends in a leap.
        [d1, d2] = iauCal2jd(iy, im, id);
        as1 = d1;
        as2 = d2 - 1 + datd;
        
        % Is the TAI after this point?
        da = a1 - as1;
        da = da + ( a2 - as2 );
        if ( da > 0 )
            % Yes:  fraction of the current UTC day that has elapsed.
            fd = da * DAYSEC / ( DAYSEC + ddats );
            
            % Ramp TAI-UTC to bring about SOFA's JD(UTC) convention.
            if(fd <= 1)
                temp = fd;
            else
                temp = 1;
            end
            datd = datd + ddats * temp / DAYSEC;
        end
        
        % Done.
        break
    end
    dats1 = dats2;
end

% Subtract the (possibly adjusted) TAI-UTC from TAI to give UTC.
a2 = a2 - datd;

% Return the UTC result, preserving the TAI order.
if (big1)
    utc1 = a1;
    utc2 = a2;
else
    utc1 = a2;
    utc2 = a1;
end

