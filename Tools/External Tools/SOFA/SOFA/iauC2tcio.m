%  - - - - - - - - - -
%   i a u C 2 t c i o
%  - - - - - - - - - -
%
%  Assemble the celestial to terrestrial matrix from CIO-based
%  components (the celestial-to-intermediate matrix, the Earth Rotation
%  Angle and the polar motion matrix).
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     rc2i               celestial-to-intermediate matrix
%     era                Earth rotation angle
%     rpom               polar-motion matrix
%
%  Returned:
%     rc2t               celestial-to-terrestrial matrix
%
%  Notes:
%  1) This function constructs the rotation matrix that transforms
%     vectors in the celestial system into vectors in the terrestrial
%     system.  It does so starting from precomputed components, namely
%     the matrix which rotates from celestial coordinates to the
%     intermediate frame, the Earth rotation angle and the polar motion
%     matrix.  One use of the present function is when generating a
%     series of celestial-to-terrestrial matrices where only the Earth
%     Rotation Angle changes, avoiding the considerable overhead of
%     recomputing the precession-nutation more often than necessary to
%     achieve given accuracy objectives.
%
%  2) The relationship between the arguments is as follows:
%
%        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
%
%              = rc2t * [CRS]
%
%     where [CRS] is a vector in the Geocentric Celestial Reference
%     System and [TRS] is a vector in the International Terrestrial
%     Reference System (see IERS Conventions 2003).
%
%  Called:
%     iauCr        copy r-matrix
%     iauRz        rotate around Z-axis
%     iauRxr       product of two r-matrices
%
%  Reference:
%     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
%     IERS Technical Note No. 32, BKG
%
%  This revision:  2008 May 11
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rc2t = iauC2tcio(rc2i, era, rpom)

% Construct the matrix.
r = iauCr(rc2i);
r = iauRz(era, r);
rc2t = iauRxr(rpom, r);

