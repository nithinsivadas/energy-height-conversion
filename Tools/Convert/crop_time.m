function [ yNew,timeNew ] = crop_time( y,time,timeMin, timeMax )
%% crop_time.m Crop a 2D matrix NxM, with time dimension M to a new time dimension M*<M
%--------------------------------------------------------------------------
% Input
%------
%  y       : 2D matrix [NxM], where M is the time dimension
%  time    : Time Array 
%  timeMin :  Initial crop time in matlab time units
%  timeMax :  Final crop time in matlab time units
%--------------------------------------------------------------------------
% Outut
%------
%  yNew    : Cropped 2D matrix [NxM], where M* is the cropped time dimension
%  timeNew : Cropped time Array 

%--------------------------------------------------------------------------
% Modified: 21st Sep 2016| 25th Jan 2017 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%--------------------------------------------------------------------------
%%
	nTime = length(time);
%     yTemp = reshape(y,nTime,[]); % not sure if this works as I want it to. 
    yTemp = y;
	timeTemp = time;

	if nargin < 4
	    timeMax = max(timeTemp);
	    itimeMax = length(timeTemp);
	else
		itimeMax=find_time(timeTemp, timeMax);
	end

	if nargin < 3
	    timeMin = min(timeTemp);
	    itimeMin = 1;
	else
	    itimeMin = find_time(timeTemp, timeMin);

    end
    
    if length(timeTemp)>1
        % Cropping the time array and density matrix to that prescribed by the user
        timeNew = timeTemp(itimeMin:itimeMax); % Changed itimeMin+1 -> itimeMin on 25th Jan 2017
        yNew = yTemp(itimeMin:itimeMax,:);
    else
        timeNew = timeTemp;
        yNew = yTemp;
    end
    
% 	[isThereNAN, totalNAN] = check_nan(yNew);

end
