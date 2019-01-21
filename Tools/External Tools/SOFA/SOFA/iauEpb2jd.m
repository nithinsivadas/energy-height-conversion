%  - - - - - - - - - -
%   i a u E p b 2 j d
%  - - - - - - - - - -
%
%  Besselian Epoch to Julian Date.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     epb      Besselian Epoch (e.g. 1957.3D0)
%
%  Returned:
%     djm0     MJD zero-point: always 2400000.5
%     djm      Modified Julian Date
%
%  Note:
%     The Julian Date is returned in two pieces, in the usual SOFA
%     manner, which is designed to preserve time resolution.  The
%     Julian Date is available as a single number by adding djm0 and
%     djm.
%
%  Reference:
%     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
%
%  This revision:  2008 May 24
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [djm0, djm] = iauEpb2jd(epb)

constants

djm0 = 2400000.5;
djm  = 15019.81352 + (epb - 1900) * DTY;

