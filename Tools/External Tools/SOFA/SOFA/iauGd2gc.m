%  - - - - - - - - -
%   i a u G d 2 g c
%  - - - - - - - - -
%
%  Transform geodetic coordinates to geocentric using the specified
%  reference ellipsoid.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical transformation.
%
%  Given:
%     n            ellipsoid identifier (Note 1)
%     elong        longitude (radians, east +ve)
%     phi          latitude (geodetic, radians, Note 3)
%     height       height above ellipsoid (geodetic, Notes 2,3)
%
%  Returned:
%     xyz          geocentric vector (Note 2)
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
%  2) The height (height, given) and the geocentric vector (xyz,
%     returned) are in meters.
%
%  3) No validation is performed on the arguments elong, phi and
%     height.  An error status -1 means that the identifier n is
%     illegal.  An error status -2 protects against cases that would
%     lead to arithmetic exceptions.  In all error cases, xyz is set
%     to zeros.
%
%  4) The inverse transformation is performed in the function iauGc2gd.
%
%  Called:
%     iauEform     Earth reference ellipsoids
%     iauGd2gce    geodetic to geocentric transformation, general
%     iauZp        zero p-vector
%
%  This revision:  2012 February 23
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xyz = iauGd2gc(n, elong, phi, height)

% Obtain reference ellipsoid parameters.
[a, f] = iauEform(n);

% transform longitude, geodetic latitude, height to x,y,z.
xyz = iauGd2gce(a, f, elong, phi, height);

