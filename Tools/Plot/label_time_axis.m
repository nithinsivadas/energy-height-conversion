function [TTick,TTickLim]=label_time_axis(time, setLabel, timeTick, timeMinStr, timeMaxStr )
%% LABEL_TIME_AXIS Sets the time axis (x-axis) time-tick locations and labels
%--------------------------------------------------------------------------
% Input:
%-------
% time          : 1xN - A time array [matlab units]
% setLabel      : true or false, if true there will be  tick marks, and labels
%                 if false there will be only tick marks
% timeTick      : Value of time interval between tick marks [Hr]
% timeMinStr    : Minimum value of time [String] Ex: '26 Mar 2008 11:00'
% timeMaxStr    : Maximum value of time [String] Ex: '26 Mar 2008 11:00'

%----------------------------------------------------------------------------
% Modified: 21st Sep 2016
%         : 29th Sep 2016 : Included label_time_axis function to shorten
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

%% Placing and calculating time tick marks
%% Default Settings
	if nargin>5
	    error('plot_eflux:TooManyInputs','required at most 5 inputs');
    end;
    if nargin<5
        timeMaxStr=datestr(max(time));
    end;
    if nargin<4
       timeMinStr=datestr(min(time));
    end;
    if nargin<3
	    timeTick=0.5;
    end;
    if nargin<2
	    setLabel=true;
    end;
    if nargin<1
	    error('plot_eflux:TooFewInputs','required at least 1 inputs');
	end;
    
    [TTick, DateNumBeg, DateNumEnd] = get_axes_time_tick_values( time, timeTick, timeMinStr, timeMaxStr );
%     datestr(DateNumBeg)
%     datestr(DateNumEnd)
%     size(TTick)
%% Creating X-axis Ticks & Labels
 	set(gca,'XTick',TTick);
    
    if setLabel == true
        set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
        delta = -0.11;
        txtHandle = text(-0.06, delta, 'UT [HH:MM]','Units','normalized','FontSize',8);
        set(txtHandle, 'HorizontalAlignment','right','VerticalAlignment','middle');
    elseif setLabel == false
        set(gca,'XTickLabel','');
    end;
    
    TTickLim=[DateNumBeg DateNumEnd];
    xlim(TTickLim);
    
end

