%  - - - - - - - - - -
%   i a u E e c t 0 0
%  - - - - - - - - - -
%
%  Equation of the equinoxes complementary terms, consistent with
%  IAU 2000 resolutions.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     date1,date2     TT as a 2-part Julian Date (Note 1)
%
%  Returned (function value):
%                     complementary terms (Note 2)
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
%  2) The "complementary terms" are part of the equation of the
%     equinoxes (EE), classically the difference between apparent and
%     mean Sidereal Time:
%
%        GAST = GMST + EE
%
%     with:
%
%        EE = dpsi * cos(eps)
%
%     where dpsi is the nutation in longitude and eps is the obliquity
%     of date.  However, if the rotation of the Earth were constant in
%     an inertial frame the classical formulation would lead to
%     apparent irregularities in the UT1 timescale traceable to side-
%     effects of precession-nutation.  In order to eliminate these
%     effects from UT1, "complementary terms" were introduced in 1994
%     (IAU, 1994) and took effect from 1997 (Capitaine and Gontier,
%     1993):
%
%        GAST = GMST + CT + EE
%
%     By convention, the complementary terms are included as part of
%     the equation of the equinoxes rather than as part of the mean
%     Sidereal Time.  This slightly compromises the "geometrical"
%     interpretation of mean sidereal time but is otherwise
%     inconsequential.
%
%     The present function computes CT in the above expression,
%     compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
%     IERS Conventions 2003).
%
%  Called:
%     iauFal03     mean anomaly of the Moon
%     iauFalp03    mean anomaly of the Sun
%     iauFaf03     mean argument of the latitude of the Moon
%     iauFad03     mean elongation of the Moon from the Sun
%     iauFaom03    mean longitude of the Moon's ascending node
%     iauFave03    mean longitude of Venus
%     iauFae03     mean longitude of Earth
%     iauFapa03    general accumulated precession in longitude
%
%  References:
%     Capitaine, N. & Gontier, A.-M., Astron. Astrophys., 275,
%     645-650 (1993)
%     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
%     implement the IAU 2000 definition of UT1", Astronomy &
%     Astrophysics, 406, 1135-1149 (2003)
%     IAU Resolution C7, Recommendation 3 (1994)
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2009 December 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eect = iauEect00(date1, date2)

constants

% t:  Time since J2000.0, in Julian centuries
% fa: Fundamental arguments

% -----------------------------------------
% The series for the EE complementary terms
% -----------------------------------------

% nfa:  coefficients of l,l',F,D,Om,LVe,LE,pA
% s, c: sine and cosine coefficients

% Terms of order t^0
e0 = [
   % 1-10
      0,  0,  0,  0,  1,  0,  0,  0, 2640.96e-6, -0.39e-6 ;
      0,  0,  0,  0,  2,  0,  0,  0,   63.52e-6, -0.02e-6 ;
      0,  0,  2, -2,  3,  0,  0,  0,   11.75e-6,  0.01e-6 ;
      0,  0,  2, -2,  1,  0,  0,  0,   11.21e-6,  0.01e-6 ;
      0,  0,  2, -2,  2,  0,  0,  0,   -4.55e-6,  0.00e-6 ;
      0,  0,  2,  0,  3,  0,  0,  0,    2.02e-6,  0.00e-6 ;
      0,  0,  2,  0,  1,  0,  0,  0,    1.98e-6,  0.00e-6 ;
      0,  0,  0,  0,  3,  0,  0,  0,   -1.72e-6,  0.00e-6 ;
      0,  1,  0,  0,  1,  0,  0,  0,   -1.41e-6, -0.01e-6 ;
      0,  1,  0,  0, -1,  0,  0,  0,   -1.26e-6, -0.01e-6 ;
   % 11-20
      1,  0,  0,  0, -1,  0,  0,  0,   -0.63e-6,  0.00e-6 ;
      1,  0,  0,  0,  1,  0,  0,  0,   -0.63e-6,  0.00e-6 ;
      0,  1,  2, -2,  3,  0,  0,  0,    0.46e-6,  0.00e-6 ;
      0,  1,  2, -2,  1,  0,  0,  0,    0.45e-6,  0.00e-6 ;
      0,  0,  4, -4,  4,  0,  0,  0,    0.36e-6,  0.00e-6 ;
      0,  0,  1, -1,  1, -8, 12,  0,   -0.24e-6, -0.12e-6 ;
      0,  0,  2,  0,  0,  0,  0,  0,    0.32e-6,  0.00e-6 ;
      0,  0,  2,  0,  2,  0,  0,  0,    0.28e-6,  0.00e-6 ;
      1,  0,  2,  0,  3,  0,  0,  0,    0.27e-6,  0.00e-6 ;
      1,  0,  2,  0,  1,  0,  0,  0,    0.26e-6,  0.00e-6 ;
   % 21-30
      0,  0,  2, -2,  0,  0,  0,  0,   -0.21e-6,  0.00e-6 ;
      0,  1, -2,  2, -3,  0,  0,  0,    0.19e-6,  0.00e-6 ;
      0,  1, -2,  2, -1,  0,  0,  0,    0.18e-6,  0.00e-6 ;
      0,  0,  0,  0,  0,  8,-13, -1,   -0.10e-6,  0.05e-6 ;
      0,  0,  0,  2,  0,  0,  0,  0,    0.15e-6,  0.00e-6 ;
      2,  0, -2,  0, -1,  0,  0,  0,   -0.14e-6,  0.00e-6 ;
      1,  0,  0, -2,  1,  0,  0,  0,    0.14e-6,  0.00e-6 ;
      0,  1,  2, -2,  2,  0,  0,  0,   -0.14e-6,  0.00e-6 ;
      1,  0,  0, -2, -1,  0,  0,  0,    0.14e-6,  0.00e-6 ;
      0,  0,  4, -2,  4,  0,  0,  0,    0.13e-6,  0.00e-6 ;
   % 31-33
      0,  0,  2, -2,  4,  0,  0,  0,   -0.11e-6,  0.00e-6 ;
      1,  0, -2,  0, -3,  0,  0,  0,    0.11e-6,  0.00e-6 ;
      1,  0, -2,  0, -1,  0,  0,  0,    0.11e-6,  0.00e-6
   ];

% Terms of order t^1
e1 = [0,  0,  0,  0,  1,  0,  0,  0,    -0.87e-6,  0.00e-6];

% Number of terms in the series
[NE0, temp] = size(e0);
[NE1, temp] = size(e1);

%--------------------------------------------------------------------

% Interval between fundamental epoch J2000.0 and current date (JC).
t = ((date1 - DJ00) + date2)/DJC;

% Fundamental Arguments (from IERS Conventions 2003)

% Mean anomaly of the Moon.
fa(1) = iauFal03(t);

% Mean anomaly of the Sun.
fa(2) = iauFalp03(t);

% Mean longitude of the Moon minus that of the ascending node.
fa(3) = iauFaf03(t);

% Mean elongation of the Moon from the Sun.
fa(4) = iauFad03(t);

% Mean longitude of the ascending node of the Moon.
fa(5) = iauFaom03(t);

% Mean longitude of Venus.
fa(6) = iauFave03(t);

% Mean longitude of Earth.
fa(7) = iauFae03(t);

% General precession in longitude.
fa(8) = iauFapa03(t);

% Evaluate the EE complementary terms.
s0 = 0;
s1 = 0;

for i = NE0:-1:1
    a = 0;
    for j = 1:8
        a = a + e0(i,j) * fa(j);
    end
    s0 = s0 + e0(i,9) * sin(a) + e0(i,10) * cos(a);
end

for i = NE1:-1:1
    a = 0;
    for j = 1:8
        a = a + e1(i,j) * fa(j);
    end
    s1 = s1 + e1(i,9) * sin(a) + e1(i,10) * cos(a);
end

eect = (s0 + s1 * t ) * DAS2R;

