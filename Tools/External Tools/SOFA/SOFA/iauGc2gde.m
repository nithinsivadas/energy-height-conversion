%  - - - - - - - - - -
%   i a u G c 2 g d e
%  - - - - - - - - - -
%
%  Transform geocentric coordinates to geodetic for a reference
%  ellipsoid of specified form.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     a            equatorial radius (Notes 2,4)
%     f            flattening (Note 3)
%     xyz(3)       geocentric vector (Note 4)
%
%  Returned:
%     elong        longitude (radians, east +ve)
%     phi          latitude (geodetic, radians)
%     height       height above ellipsoid (geodetic, Note 4)
%
%  Notes:
%  1) This function is based on the GCONV2H Fortran subroutine by
%     Toshio Fukushima (see reference).
%
%  2) The equatorial radius, a, can be in any units, but meters is
%     the conventional choice.
%
%  3) The flattening, f, is (for the Earth) a value around 0.00335,
%     i.e. around 1/298.
%
%  4) The equatorial radius, a, and the geocentric vector, xyz,
%     must be given in the same units, and determine the units of
%     the returned height, height.
%
%  5) If an error occurs (status < 0), elong, phi and height are
%     unchanged.
%
%  6) The inverse transformation is performed in the function
%     iauGd2gce.
%
%  7) The transformation for a standard ellipsoid (such as WGS84) can
%     more conveniently be performed by calling iauGc2gd, which uses a
%     numerical code to identify the required A and F values.
%
%  Reference:
%     Fukushima, T., "Transformation from Cartesian to geodetic
%     coordinates accelerated by Halley's method", J.Geodesy (2006)
%     79: 689-693
%
%  This revision:  2012 February 23
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elong, phi, height] = iauGc2gde(a, f, xyz)

constants

% Functions of ellipsoid parameters (with further validation of f).
aeps2 = a*a * 1e-32;
e2 = (2 - f) * f;
e4t = e2*e2 * 1.5;
ec2 = 1 - e2;
ec = sqrt(ec2);
b = a * ec;

% Cartesian components.
x = xyz(1);
y = xyz(2);
z = xyz(3);

% Distance from polar axis squared.
p2 = x*x + y*y;

% Longitude.
if (p2 ~= 0)
    elong = atan2(y, x);
else
    elong = 0;
end

% Unsigned z-coordinate.
absz = abs(z);

% Proceed unless polar case.
if (p2 > aeps2)
    % Distance from polar axis.
    p = sqrt(p2);
    
    % Normalization.
    s0 = absz / a;
    pn = p / a;
    zc = ec * s0;
    
    % Prepare Newton correction factors.
    c0 = ec * pn;
    c02 = c0 * c0;
    c03 = c02 * c0;
    s02 = s0 * s0;
    s03 = s02 * s0;
    a02 = c02 + s02;
    a0 = sqrt(a02);
    a03 = a02 * a0;
    d0 = zc*a03 + e2*s03;
    f0 = pn*a03 - e2*c03;
    
    % Prepare Halley correction factor.
    b0 = e4t * s02 * c02 * pn * (a0 - ec);
    s1 = d0*f0 - b0*s0;
    cc = ec * (f0*f0 - b0*c0);
    
    % Evaluate latitude and height.
    phi = atan(s1/cc);
    s12 = s1 * s1;
    cc2 = cc * cc;
    height = (p*cc + absz*s1 - a * sqrt(ec2*s12 + cc2))/sqrt(s12 + cc2);
else    
    % Exception: pole.
    phi = DPI / 2;
    height = absz - b;
end

% Restore sign of latitude.
if (z < 0)
    phi = -phi;
end

