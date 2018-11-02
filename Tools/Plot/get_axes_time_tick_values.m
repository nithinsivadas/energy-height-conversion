function [ TTick, DateNumBeg, DateNumEnd ] = get_axes_time_tick_values( time, timeTick, timeMinStr, timeMaxStr )
%% get_axes_time_tick_values.m Calculates the time tick values in matlab units 
%--------------------------------------------------------------------------
% Input:
%--------
% time     : Time array in matlab units [1xN]
% timeTick : Time tick interval in hr: example - 0.5 will produce half hour
%            intervals
% timeMinStr: The beginning time value [in string Example: 26 Mar 2008
%             11:00]
% timeMaxStr: The end time value [in string]
%--------------------------------------------------------------------------
% Output:
%--------
% TTick     : Matlab time values where a axis tick ought to be displayed
% DateNumBeg: Lower time limit of the axes
% DateNumEnd: Upper time limit of the axes
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------

    itimeStart = find_time(time, timeMinStr);
    DateNumBeg = time(itimeStart);
    
    itimeEnd = find_time(time, timeMaxStr);
    DateNumEnd = time(itimeEnd);

	dt=timeTick;
	TTick=(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):...
        dt/24:...
        floor(DateNumEnd)+(ceil((DateNumEnd-floor(DateNumEnd))/(dt/24))/(24/dt)));
    
     if TTick(1) > DateNumBeg+0.5*dt/24
         TTick(2:end+1)=TTick;
         TTick(1)=DateNumBeg;
     else
         TTick(1)=DateNumBeg;
     end
     
     if TTick(end)<DateNumEnd-0.5*dt/24 && TTick(end)<DateNumEnd
        TTick(end+1)=DateNumEnd;
     elseif DateNumEnd>TTick(end-1)+0.5*dt/24
        TTick(end)=DateNumEnd;
     else
         TTick(end)=[];
     end
     
end

