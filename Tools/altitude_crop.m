function [yNew, altNew] = altitude_crop(Y, alt, minAlt, maxAlt)
%%altitudeCrop.m Crop Matrix Y with N time values and M altitude values to M* altitude values

% Input
%  y     : 2D matrix [NxM], where N is the time dimension and is monotonically
%         increasing; while M is the altitude dimension
%  alt    : Altitude array
%  minAlt : Lower limit of cropped altitude dimension 
%  maxAlt : Upper limit of cropped altitude dimension

% Output
% yNew    : 2D matrix [NxM*], where N is the time dimension and M* is the cropped altitude
% 		    dimension
% altNew  : The cropped altitude array
%
%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

% Storing orginal values in a temporary array
	yTemp   = y;
	altTemp = alt;

% Assigning default values
    switch nargin
    	case 4
    		[c,AltEndIndex] = min(abs(altTemp-maxAlt));
    		[c,AltBegIndex] = min(abs(altTemp-minAlt));		
    	case 3
        	maxAlt=max(altTemp);
    		AltEndIndex=length(altTemp);
    		[c,AltBegIndex] = min(abs(altTemp-minAlt));
    	case 2
    		maxAlt=max(altTemp);
    		AltEndIndex=length(altTemp);
    		minAlt=min(altTemp);
    		AltBegIndex = 1;	
    	otherwise
    		if nargin>4
    			error('TooManyInputs: Number of input variables should be at most 4')
    		else
    			error('TooFewInputs: Number of input variables should at least be 2');
    		end
    end

% Cropping the altitude array and density matrix to that prescribed by the user
altNew	= altTemp(AltBegIndex+1:AltEndIndex);
yNew    = yTemp(AltBegIndex+1:AltEndIndex,:);

end