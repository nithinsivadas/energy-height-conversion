%  - - - - - - -
%   i a u E p b
%  - - - - - - -
%
%  Julian Date to Besselian Epoch.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     dj1,dj2    Julian Date (see note)
%
%  Returned (function value):
%                Besselian Epoch.
%
%  Note:
%     The Julian Date is supplied in two pieces, in the usual SOFA
%     manner, which is designed to preserve time resolution.  The
%     Julian Date is available as a single number by adding dj1 and
%     dj2.  The maximum resolution is achieved if dj1 is 2451545D0
%     (J2000.0).
%
%  Reference:
%     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
%
%  This revision:  2009 December 16
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epb = iauEpb(dj1, dj2)

constants

% J2000.0 minus B1900.0 (2415019.81352) in Julian days
D1900 = 36524.68648;

epb = 1900 + ((dj1 - DJ00) + (dj2 + D1900)) / DTY;

