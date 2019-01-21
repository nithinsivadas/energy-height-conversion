%  - - - - - - - - -
%   i a u F w 2 x y
%  - - - - - - - - -
%
%  CIP X,Y given Fukushima-Williams bias-precession-nutation angles.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     gamb         F-W angle gamma_bar (radians)
%     phib         F-W angle phi_bar (radians)
%     psi          F-W angle psi (radians)
%     eps          F-W angle epsilon (radians)
%
%  Returned:
%     x,y          CIP X,Y ("radians")
%
%  Notes:
%  1) Naming the following points:
%
%           e = J2000.0 ecliptic pole,
%           p = GCRS pole
%           E = ecliptic pole of date,
%     and   P = CIP,
%
%     the four Fukushima-Williams angles are as follows:
%
%        gamb = gamma = epE
%        phib = phi = pE
%        psi = psi = pEP
%        eps = epsilon = EP
%
%  2) The matrix representing the combined effects of frame bias,
%     precession and nutation is:
%
%        NxPxB = R_1(-epsA).R_3(-psi).R_1(phib).R_3(gamb)
%
%     X,Y are elements (3,1) and (3,2) of the matrix.
%
%  Called:
%     iauFw2m      F-W angles to r-matrix
%     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
%
%  Reference:
%     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
%
%  This revision:  2009 December 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = iauFw2xy(gamb, phib, psi, eps)

% Form NxPxB matrix.
r = iauFw2m(gamb, phib, psi, eps);

% Extract CIP X,Y.
[x, y] = iauBpn2xy(r);

