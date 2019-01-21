%  - - - - - - - -
%   i a u B p 0 0
%  - - - - - - - -
%
%  Frame bias and precession, IAU 2000.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  canonical model.
%
%  Given:
%     date1,date2           TT as a 2-part Julian Date (Note 1)
%
%  Returned:
%     rb                    frame bias matrix (Note 2)
%     rp                    precession matrix (Note 3)
%     rbp                   bias-precession matrix (Note 4)
%
%  Notes:
%
%  1) The TT date date1+date2 is a Julian Date, apportioned in any
%     convenient way between the two arguments.  For example,
%     JD(TT)=2450123.7 could be expressed in any of these ways,
%     among others:
%
%             date1         date2
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
%  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
%     applying frame bias.
%
%  3) The matrix rp transforms vectors from J2000.0 mean equator and
%     equinox to mean equator and equinox of date by applying
%     precession.
%
%  4) The matrix rbp transforms vectors from GCRS to mean equator and
%     equinox of date by applying frame bias then precession.  It is
%     the product rp x rb.
%
%  5) It is permissible to re-use the same array in the returned
%     arguments.  The arrays are filled in the order given.
%
%  Called:
%     iauBi00      frame bias components, IAU 2000
%     iauPr00      IAU 2000 precession adjustments
%     iauIr        initialize r-matrix to identity
%     iauRx        rotate around X-axis
%     iauRy        rotate around Y-axis
%     iauRz        rotate around Z-axis
%     iauCr        copy r-matrix
%     iauRxr       product of two r-matrices
%
%  Reference:
%     "Expressions for the Celestial Intermediate Pole and Celestial
%     Ephemeris Origin consistent with the IAU 2000A precession-
%     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
%
%     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
%          intermediate origin" (CIO) by IAU 2006 Resolution 2.
%
%  This revision:  2010 January 18
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rb, rp, rbp] = iauBp00(date1, date2)

constants

% J2000.0 obliquity (Lieske et al. 1977)
EPS0 = 84381.448 * DAS2R;

% Interval between fundamental epoch J2000.0 and current date (JC).
t = ((date1 - DJ00) + date2) / DJC;

% Frame bias.
[dpsibi, depsbi, dra0] = iauBi00;

% Precession angles (Lieske et al. 1977)
psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * DAS2R;
oma77  =       EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * DAS2R;
chia   = (  10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * DAS2R;

% Apply IAU 2000 precession corrections.
[dpsipr, depspr] = iauPr00(date1, date2);
psia = psia77 + dpsipr;
oma  = oma77  + depspr;

% Frame bias matrix: GCRS to J2000.0.
rbw = zeros(3);
rbw = iauIr(rbw);
rbw = iauRz(dra0, rbw);
rbw = iauRy(dpsibi * sin(EPS0), rbw);
rbw = iauRx(-depsbi, rbw);
rb = iauCr(rbw);

% Precession matrix: J2000.0 to mean of date.
rp = zeros(3);
rp = iauIr(rp);
rp = iauRx(EPS0, rp);
rp = iauRz(-psia, rp);
rp = iauRx(-oma, rp);
rp = iauRz(chia, rp);

% Bias-precession matrix: GCRS to mean of date.
rbp = iauRxr(rp, rbw);

