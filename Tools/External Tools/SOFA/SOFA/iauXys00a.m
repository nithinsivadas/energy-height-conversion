%  - - - - - - - - - -
%   i a u X y s 0 0 a
%  - - - - - - - - - -
%
%  For a given TT date, compute the X,Y coordinates of the Celestial
%  Intermediate Pole and the CIO locator s, using the IAU 2000A
%  precession-nutation model.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     date1,date2     TT as a 2-part Julian Date (Note 1)
%
%  Returned:
%     x,y             Celestial Intermediate Pole (Note 2)
%     s               the CIO locator s (Note 2)
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
%  2) The Celestial Intermediate Pole coordinates are the x,y
%     components of the unit vector in the Geocentric Celestial
%     Reference System.
%
%  3) The CIO locator s (in radians) positions the Celestial
%     Intermediate Origin on the equator of the CIP.
%
%  4) A faster, but slightly less accurate result (about 1 mas for
%     X,Y), can be obtained by using instead the iauXys00b function.
%
%  Called:
%     iauPnm00a    classical NPB matrix, IAU 2000A
%     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
%     iauS00       the CIO locator s, given X,Y, IAU 2000A
%
%  Reference:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2008 May 12
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, s] = iauXys00a(date1, date2)

% Form the bias-precession-nutation matrix, IAU 2000A.
rbpn = iauPnm00a(date1, date2);

% Extract X,Y.
[x y] = iauBpn2xy(rbpn);

% Obtain s.
s = iauS00(date1, date2, x, y);

