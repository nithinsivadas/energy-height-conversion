%  - - - - - - - - -
%   i a u E e 0 0 a
%  - - - - - - - - -
%
%  Equation of the equinoxes, compatible with IAU 2000 resolutions.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     date1,date2      TT as a 2-part Julian Date (Note 1)
%
%  Returned (function value):
%                      equation of the equinoxes (Note 2)
%
%  Notes:
%  1) The TT date date1+date2 is a Julian Date, apportioned in any
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
%  3) The result is compatible with the IAU 2000 resolutions.  For
%     further details, see IERS Conventions 2003 and Capitaine et al.
%     (2002).
%
%  Called:
%     iauPr00      IAU 2000 precession adjustments
%     iauObl80     mean obliquity, IAU 1980
%     iauNut00a    nutation, IAU 2000A
%     iauEe00      equation of the equinoxes, IAU 2000
%
%  References:
%     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
%     implement the IAU 2000 definition of UT1", Astronomy &
%     Astrophysics, 406, 1135-1149 (2003).
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004).
%
%  This revision:  2008 May 16
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ee = iauEe00a(date1, date2)

% IAU 2000 precession-rate adjustments.
[dpsipr, depspr] = iauPr00(date1, date2);

% Mean obliquity, consistent with IAU 2000 precession-nutation.
epsa = iauObl80(date1, date2) + depspr;

% Nutation in longitude.
[dpsi, deps] = iauNut00a(date1, date2);

% Equation of the equinoxes.
ee = iauEe00(date1, date2, epsa, dpsi);

