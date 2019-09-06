function [ rotatedArray ] = rotate_array( array, alpha )
%% rotate_array Rotate the array clockwise w.r.t data (and counter-clockwise w.r.t to previous configuration)
%--------------------------------------------------------------------------
%Input
%-----
% array - any array consisting of angles in degrees [nArray]
% alpha - the angle through which the elements of the array are to be rotated [in deg]
%--------------------------------------------------------------------------
% Output
%-------
% rotatedArray - rotated array in degrees [nArray]
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------

a1 = array - alpha;
% I = find(a1<0);
% a1(I) = 360 + a1(I);
rotatedArray = wrapTo360(a1);
end

