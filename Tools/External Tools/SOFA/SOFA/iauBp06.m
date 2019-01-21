%  - - - - - - - -
%   i a u B p 0 6
%  - - - - - - - -
%
%  Frame bias and precession, IAU 2006.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
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
%  3) The matrix rp transforms vectors from mean J2000.0 to mean of
%     date by applying precession.
%
%  4) The matrix rbp transforms vectors from GCRS to mean of date by
%     applying frame bias then precession.  It is the product rp x rb.
%
%  Called:
%     iauPfw06     bias-precession F-W angles, IAU 2006
%     iauFw2m      F-W angles to r-matrix
%     iauPmat06    PB matrix, IAU 2006
%     iauTr        transpose r-matrix
%     iauRxr       product of two r-matrices
%
%  References:
%     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
%     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
%
%  This revision:  2009 December 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rb, rp, rbp] = iauBp06(date1, date2)

% B matrix.
[gamb, phib, psib, epsa] = iauPfw06(DJM0, DJM00);
rb = iauFw2m(gamb, phib, psib, epsa);

% PxB matrix.
rbp = iauPmat06(date1, date2);

% P matrix.
rbt = iauTr(rb);
rp = iauRxr(rbp, rbt);

