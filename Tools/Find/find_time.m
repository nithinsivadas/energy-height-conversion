function [timeNo] = find_time( time, thisTimeStr)

%%find_time.m Find the index of a 1-D time array, given 
%           a particular time value 
% Input
%  time         : 1-D time Array 
%  thisTimeStr	: String identifying the time to be found 
%                from the array time '26 Mar 2008 11:00'
% Output:
%  timeNo       : The index of the time array that points 
%           	 to the time specified by thisTime 

%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

	[c,timeNo] = min(abs(time-datenum(thisTimeStr)));

end
