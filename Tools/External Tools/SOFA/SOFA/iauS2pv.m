%  - - - - - - - -
%   i a u S 2 p v
%  - - - - - - - -
%
%  Convert position/velocity from spherical to Cartesian coordinates.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     theta              longitude angle (radians)
%     phi                latitude angle (radians)
%     r                  radial distance
%     td                 rate of change of theta
%     pd                 rate of change of phi
%     rd                 rate of change of r
%
%  Returned:
%     pv(2,3)            pv-vector
%
%  This revision:  2008 May 25
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pv = iauS2pv(theta, phi, r, td, pd, rd)

st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);
rcp = r * cp;
x = rcp * ct;
y = rcp * st;
rpd = r * pd;
w = rpd*sp - cp*rd;

pv(1,1) = x;
pv(1,2) = y;
pv(1,3) = r * sp;
pv(2,1) = -y*td - w*ct;
pv(2,2) =  x*td - w*st;
pv(2,3) = rpd*cp + sp*rd;

