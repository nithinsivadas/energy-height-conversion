%  - - - - - - - - - -
%   i a u C 2 t 0 0 a
%  - - - - - - - - - -
%
%  Form the celestial to terrestrial matrix given the date, the UT1 and
%  the polar motion, using the IAU 2000A nutation model.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     tta,ttb           TT as a 2-part Julian Date (Note 1)
%     uta,utb           UT1 as a 2-part Julian Date (Note 1)
%     xp,yp             coordinates of the pole (radians, Note 2)
%
%  Returned:
%     rc2t              celestial-to-terrestrial matrix (Note 3)
%
%  Notes:
%
%  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
%     apportioned in any convenient way between the arguments uta and
%     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
%     these ways, among others:
%
%             uta            utb
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution is
%     acceptable.  The J2000 and MJD methods are good compromises
%     between resolution and convenience.  In the case of uta,utb, the
%     date & time method is best matched to the Earth rotation angle
%     algorithm used:  maximum precision is delivered when the uta
%     argument is for 0hrs UT1 on the day in question and the utb
%     argument lies in the range 0 to 1, or vice versa.
%
%  2) The arguments xp and yp are the coordinates (in radians) of the
%     Celestial Intermediate Pole with respect to the International
%     Terrestrial Reference System (see IERS Conventions 2003),
%     measured along the meridians to 0 and 90 deg west respectively.
%
%  3) The matrix rc2t transforms from celestial to terrestrial
%     coordinates:
%
%        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
%
%              = rc2t * [CRS]
%
%     where [CRS] is a vector in the Geocentric Celestial Reference
%     System and [TRS] is a vector in the International Terrestrial
%     Reference System (see IERS Conventions 2003), RC2I is the
%     celestial-to-intermediate matrix, ERA is the Earth rotation
%     angle and RPOM is the polar motion matrix.
%
%  4) A faster, but slightly less accurate result (about 1 mas), can
%     be obtained by using instead the iauC2t00b function.
%
%  Called:
%     iauC2i00a    celestial-to-intermediate matrix, IAU 2000A
%     iauEra00     Earth rotation angle, IAU 2000
%     iauSp00      the TIO locator s', IERS 2000
%     iauPom00     polar motion matrix
%     iauC2tcio    form CIO-based celestial-to-terrestrial matrix
%
%  Reference:
%     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG (2004)
%
%  This revision:  2009 April 1
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rc2t = iauC2t00a(tta, ttb, uta, utb, xp, yp)

% Form the celestial-to-intermediate matrix for this TT (IAU 2000A).
rc2i = iauC2i00a(tta, ttb);

% Predict the Earth rotation angle for this UT1.
era = iauEra00(uta, utb);

% Estimate s'.
sp = iauSp00(tta, ttb);

% Form the polar motion matrix.
rpom = iauPom00(xp, yp, sp);

% Combine to form the celestial-to-terrestrial matrix.
rc2t = iauC2tcio(rc2i, era, rpom);

