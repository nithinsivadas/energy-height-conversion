%  - - - - - - - - - -
%   i a u _ c 2 i 0 6 a
%  - - - - - - - - - -
%
%  Form the celestial-to-intermediate matrix for a given date using the
%  IAU 2006 precession and IAU 2000A nutation models.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     date1,date2        TT as a 2-part Julian Date (Note 1)
%
%  Returned:
%     rc2i               celestial-to-intermediate matrix (Note 2)
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
%  2) The matrix rc2i is the first stage in the transformation from
%     celestial to terrestrial coordinates:
%
%        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
%
%               =  RC2T * [CRS]
%
%     where [CRS] is a vector in the Geocentric Celestial Reference
%     System and [TRS] is a vector in the International Terrestrial
%     Reference System (see IERS Conventions 2003), ERA is the Earth
%     Rotation Angle and RPOM is the polar motion matrix.
%
%  Called:
%     iauPnm06a    classical NPB matrix, IAU 2006/2000A
%     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
%     iauS06       the CIO locator s, Given X,Y, IAU 2006
%     iauC2ixys    celestial-to-intermediate matrix, Given X,Y and s
%
%  References:
%     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG
%
%  This revision:  2008 May 13
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rc2i = iauC2i06a(date1, date2)

% Obtain the celestial-to-true matrix (IAU 2006/2000A).
rbpn = iauPnm06a(date1, date2);

% Extract the X,Y coordinates.
[x, y] = iauBpn2xy(rbpn);

% Obtain the CIO locator.
s = iauS06(date1, date2, x, y);

% Form the celestial-to-intermediate matrix.
rc2i = iauC2ixys(x, y, s);

