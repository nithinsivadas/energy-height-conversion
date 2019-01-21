%  - - - - - - - - -
%   i a u G c 2 g d
%  - - - - - - - - -
%
%  Transform geocentric coordinates to geodetic using the specified
%  reference ellipsoid.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical transformation.
%
%  Given:
%     n            ellipsoid identifier (Note 1)
%     xyz(3)       geocentric vector (Note 2)
%
%  Returned:
%     elong        longitude (radians, east +ve)
%     phi          latitude (geodetic, radians, Note 3)
%     height       height above ellipsoid (geodetic, Notes 2,3)
%
%  Notes:
%  1) The identifier n is a number that specifies the choice of
%     reference ellipsoid.  The following are supported:
%
%        n    ellipsoid
%
%        1     WGS84
%        2     GRS80
%        3     WGS72
%
%     The n value has no significance outside the SOFA software.  For
%     convenience, symbols WGS84 etc. are defined in sofam.h.
%
%  2) The geocentric vector (xyz, given) and height (height, returned)
%     are in meters.
%
%  3) An error status -1 means that the identifier n is illegal.  An
%     error status -2 is theoretically impossible.  In all error cases,
%     phi and height are both set to -1e9.
%
%  4) The inverse transformation is performed in the function iauGd2gc.
%
%  Called:
%     iauEform     Earth reference ellipsoids
%     iauGc2gde    geocentric to geodetic transformation, general
%
%  This revision:  2012 Febuary 23
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elong, phi, height] = iauGc2gd(n, xyz)

% Obtain reference ellipsoid parameters.
[a, f] = iauEform(n);

% transform x,y,z to longitude, geodetic latitude, height.
[elong, phi, height] = iauGc2gde(a, f, xyz);

