%   - - - - - - - -
%    i a u A n p m
%   - - - - - - - -
% 
%   Normalize angle into the range -pi <= a < +pi.
% 
%   This function is part of the International Astronomical Union's
%   SOFA (Standards Of Fundamental Astronomy) software collection.
% 
%   Status:  vector/matrix support function.
% 
%   Given:
%      a             angle (radians)
% 
%   Returned (function value):
%                    angle in range +/-pi
% 
%   This revision:  2008 May 16
% 
%   SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = iauAnpm(a)

constants

w = mod(a, D2PI);

if (abs(w) >= DPI)
    w = w-dsign(D2PI, a);
end

