%  - - - - - - - - - -
%   i a u E q e q 9 4
%  - - - - - - - - - -
%
%  Equation of the equinoxes, IAU 1994 model.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     date1,date2        TDB date (Note 1)
%
%  Returned (function value):
%                        equation of the equinoxes (Note 2)
%
%  Notes:
%  1) The date date1+date2 is a Julian Date, apportioned in any
%     convenient way between the two arguments.  For example,
%     JD(TT)=2450123.7 could be expressed in any of these ways,
%     among others:
%
%            date1          date2
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable.  The J2000 method is best matched to the way
%     the argument is handled internally and will deliver the
%     optimum resolution.  The MJD method and the date & time methods
%     are both good compromises between resolution and convenience.
%
%  2) The result, which is in radians, operates in the following sense:
%
%        Greenwich apparent ST = GMST + equation of the equinoxes
%
%  Called:
%     iauNut80     nutation, IAU 1980
%     iauObl80     mean obliquity, IAU 1980
%
%  References:
%     IAU Resolution C7, Recommendation 3 (1994).
%     Capitaine, N. & Gontier, A.-M., 1993, Astron. Astrophys., 275,
%     645-650.
%
%  This revision:  2008 May 24
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ee = iauEqeq94(date1, date2)

constants

% Interval between fundamental epoch J2000.0 and given date (JC).
t = ((date1 - DJ00) + date2) / DJC;

% Longitude of the mean ascending node of the lunar orbit on the
% ecliptic, measured from the mean equinox of date.
om = iauAnpm((450160.280 + (-482890.539 ...
   + (7.455 + 0.008 * t) * t) * t) * DAS2R ...
   + mod(-5.0 * t, 1) * D2PI);

% Nutation components and mean obliquity.
[dpsi, deps] = iauNut80(date1, date2);
eps0 = iauObl80(date1, date2);

% Equation of the equinoxes.
ee = dpsi*cos(eps0) + DAS2R*(0.00264*sin(om) + 0.000063*sin(om+om));

