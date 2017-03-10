function h=plot_2D_time_series(time, yAxis, zValue, timeTick, mode, timeMinStr, timeMaxStr)

%% plot2D.m Plots 2D time series data
%--------------------------------------------------------------------------
% Input:
%-------
%  time 	: 1xN - A time array [matlab units]
%  yAxis 	: Mx1 - Values of the parameter on vertical axis (y-axis)
%                 Example: Energy bin [eV]
%  zValue 	: MxN - A matrix specifying zValues along time and y-axis
%                 Example: Energy flux [eV/ cm^2 sec sr eV]
%  timeTick : Value of time interval between tick marks [Hr]
%  mode     : Uses a set of pre-determined plot options
%        '1': Assumes one is plotting Electron Density
%        '2': Assumes one is plotting Production rates
%        '3': Assumes one is plotting Differential Energy Flux
% timeMin   : Minimum value of time [String] Ex: '26 Mar 2008 11:00'
% timeMax   : Maximum value of time [String] Ex: '26 Mar 2008 11:00'
%--------------------------------------------------------------------------
% Output:
%-------
% h         : pcolor handle
%----------------------------------------------------------------------------
% Modified: 21st Sep 2016
%         : 29th Sep 2016 : Included label_time_axis function to shorten
%         : 15th Oct 2016 : Changed surf plot to pcolor
%         : 9th Mar 2017 : Modified datestr format
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

%% Default Settings
	if nargin>7
	    error('plot_eflux:TooManyInputs','required at most 7 inputs');
    end;
    if nargin<7
        timeMaxStr=datestr(max(time));
    end;
    if nargin<6
       timeMinStr=datestr(min(time));
    end;
    if nargin<5
        mode=0; % No mode
    end;
    if nargin<4
	    timeTick=0.5;
    end;
    if nargin<3
	    error('plot_eflux:TooFewInputs','required at least 3 inputs');
	end;
    
%% Creating a Surface Plot

	[X,Y]=meshgrid(time,yAxis);	

% 	surf(X,Y,zValue,'EdgeColor','none');
%   view([0 90]);    
    h=pcolor(X,Y,zValue);
    set(h,'EdgeColor','none'); 
    shading flat;
 	set(gca,'Layer','top')
	
    
    hold on;
    
%% Calculating Y axis and Color axis limits
	miny=min(yAxis); maxy=max(yAxis);
	ylim([miny maxy]);

	minz=min(min(zValue)); maxz=max(max(zValue));
	caxis([minz maxz]);
% 	colorbar;

%% Defining necessary labels according to different modes
	switch mode
		case 0
			label_time_axis(time, false, timeTick, timeMinStr, timeMaxStr);
        case 1
			set(get(gca,'YLabel'),'String','Altitude [km]');
            label_time_axis(time, true, timeTick, timeMinStr, timeMaxStr);
			title(['Electron Density log_1_0 [m^-^3] at ',datestr(time(1),'DD-mmm-YY')]);
		case 2
			set(get(gca,'YLabel'),'String','Altitude [km]');
            label_time_axis(time, true, timeTick, timeMinStr, timeMaxStr);
			title(['Production Rates log_1_0 [m^-^3 s^-^1] at ',datestr(time(1),'DD-mmm-YY')]);
		case 3
			set(get(gca,'YLabel'),'String','Energy [eV]');
            label_time_axis(time, true, timeTick, timeMinStr, timeMaxStr);
			title(['Differential Energy Flux log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1] at ',datestr(time(1),'DD-mmm-YY')]);
			set(gca,'YScale','log');
		otherwise
		    warning('No mode selected')
	end

	[isThereNAN, totalNAN] = check_nan(zValue);

end