%  - - - - - - - - -
%   i a u L d s u n
%  - - - - - - - - -
%
%  Deflection of starlight by the Sun.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     p(3)        direction from observer to star (unit vector)
%     e(3)        direction from Sun to observer (unit vector)
%     em          distance from Sun to observer (au)
%
%  Returned:
%     p1(3)       observer to deflected star (unit vector)
%
%  Notes:
%  1) The source is presumed to be sufficiently distant that its
%     directions seen from the Sun and the observer are essentially
%     the same.
%
%  2) The deflection is restrained when the angle between the star and
%     the center of the Sun is less than about 9 arcsec, falling to
%     zero for zero separation. (The chosen threshold is within the
%     solar limb for all solar-system applications.)
%
%  3) The arguments p and p1 can be the same array.
%
%  Called:
%     iauLd        light deflection by a solar-system body
%
%  This revision:   2014 September 1
%
%  SOFA release 2013-12-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p1 = iauLdsun(p, e, em)

p1 = iauLd(1, p, p, e, em, 1e-9);

