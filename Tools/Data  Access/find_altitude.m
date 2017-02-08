function [altitudeNo] = find_altitude( altitudeGrid, thisAltitude)
%% find_altitude.m Find the index of a 1-D altitude array, given a particular altitude value 

% Input
%  altitudeGrid    	: 1-D altitude Array 
%  thisaltitude     : Value identifying the index of the altitude to be found 
%                    from the array '60.5'
% Output:
%  altitudeNo: The index of the altitude array that points 
%             to the altitude specified by thisAltitude 

%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

	[c,altitudeNo] = min(abs(altitudeGrid-thisAltitude));

end