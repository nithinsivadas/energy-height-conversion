%  - - - - - - - - - -
%   i a u J d c a l f
%  - - - - - - - - - -
%
%  Julian Date to Gregorian Calendar, expressed in a form convenient
%  for formatting messages:  rounded to a specified precision.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     ndp          number of decimal places of days in fraction
%     dj1,dj2      dj1+dj2 = Julian Date (Note 1)
%
%  Returned:
%     iymdf        year, month, day, fraction in Gregorian calendar
%
%  Notes:
%
%  1) The Julian Date is apportioned in any convenient way between
%     the arguments dj1 and dj2.  For example, JD=2450123.7 could
%     be expressed in any of these ways, among others:
%
%             dj1            dj2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%  2) In early eras the conversion is from the "Proleptic Gregorian
%     Calendar";  no account is taken of the date(s) of adoption of
%     the Gregorian Calendar, nor is the AD/BC numbering convention
%     observed.
%
%  3) Refer to the function iauJd2cal.
%
%  4) NDP should be 4 or less if internal overflows are to be
%     avoided on machines which use 16-bit integers.
%
%  Called:
%     iauJd2cal    JD to Gregorian calendar
%
%  Reference:
%
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992),
%     Section 12.92 (p604).
%
%  This revision:  2010 July 27
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iymdf = iauJdcalf(ndp, dj1, dj2)

% Denominator of fraction (e.g. 100 for 2 decimal places).
if ((ndp >= 0) && (ndp <= 9))
    denom = 10.0^ndp;
else
    denom = 1.0;
end

% Copy the date, big then small, and realign to midnight.
if (dj1 >= dj2)
    d1 = dj1;
    d2 = dj2;
else
    d1 = dj2;
    d2 = dj1;
end

d2 = d2 - 0.5;

% Separate days and fractions.
f1 = mod(d1, 1.0);
f2 = mod(d2, 1.0);
d1 = floor(d1 - f1);
d2 = floor(d2 - f2);

% Round the total fraction to the specified number of places.
f = floor((f1+f2)*denom + 0.5) / denom;

% Re-assemble the rounded date and re-align to noon.
d2 = d2 + f + 0.5;

% Convert to Gregorian calendar.
[iymdf(1), iymdf(2), iymdf(3), f] = iauJd2cal(d1, d2);
iymdf(4) = int32(f * denom);

