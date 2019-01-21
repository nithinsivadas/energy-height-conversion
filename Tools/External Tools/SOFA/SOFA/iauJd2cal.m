%  - - - - - - - - - -
%   i a u J d 2 c a l
%  - - - - - - - - - -
%
%  Julian Date to Gregorian year, month, day, and fraction of a day.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     dj1,dj2   Julian Date (Notes 1, 2)
%
%  Returned (arguments):
%     iy        year
%     im        month
%     id        day
%     fd        fraction of day
%
%  Notes:
%
%  1) The earliest valid date is -68569.5 (-4900 March 1).  The
%     largest value accepted is 10^9.
%
%  2) The Julian Date is apportioned in any convenient way between
%     the arguments dj1 and dj2.  For example, JD=2450123.7 could
%     be expressed in any of these ways, among others:
%
%            dj1             dj2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%  3) In early eras the conversion is from the "proleptic Gregorian
%     calendar";  no account is taken of the date(s) of adoption of
%     the Gregorian calendar, nor is the AD/BC numbering convention
%     observed.
%
%  Reference:
%
%     Explanatory Supplement to the Astronomical Almanac,
%     P. Kenneth Seidelmann (ed), University Science Books (1992),
%     Section 12.92 (p604).
%
%  This revision:  2008 May 26
%
%  SOFA release 2012-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iy, im, id, fd] = iauJd2cal(dj1, dj2)

% ----------------- find year and days of the year ---------------
jd     = dj1 + dj2;
temp   = jd - 2415019.5;
tu     = temp / 365.25;
iy     = 1900 + floor( tu );
leapyrs= floor( ( iy-1901 )*0.25 );
days   = temp - ((iy-1900)*365.0 + leapyrs );

% ------------ check for case of beginning of a year -------------
if (days < 1.0)
    iy     = iy - 1;
    leapyrs= floor( ( iy-1901 )*0.25 );
    days   = temp - ((iy-1900)*365.0 + leapyrs );
end

% ------------------- find remaining data  -----------------------
[im,id,hr,min,sec] = days2mdh(iy, days);

if (dj1 >= dj2)
    d1 = dj1;
    d2 = dj2;
else
    d1 = dj2;
    d2 = dj1;
end

d2 = d2 - 0.5;
f1 = mod(d1, 1);
f2 = mod(d2, 1);
fd = mod(f1 + f2, 1);

if (fd < 0.0)
    fd = fd + 1;
end

