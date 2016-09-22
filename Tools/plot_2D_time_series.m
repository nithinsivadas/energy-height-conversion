function plot_2D_time_series(time, yAxis, zValue, timeTick, mode, timeMin, timeMax)

%% plot2D.m Plots the electron energy flux using an input data structure

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

%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
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
        timeMax=datestr(max(time));
    end;
    if nargin<6
       timeMin=datestr(min(time));
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

	surf(X,Y,zValue,'EdgeColor','none');
	set(gca, 'Layer','top')
	view([0 90]);

%% Placing and calculating time tick marks
    itimeStart = find_time(time, timeMin);
    DateNumBeg = time(itimeStart);
    
    itimeEnd = find_time(time, timeMax);
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
     end;
     if TTick(end)<DateNumEnd-0.5*dt/24
        TTick(end+1)=DateNumEnd;
     else
        TTick(end)=DateNumEnd;
     end;

%% Creating X-axis Ticks & Labels
 	
	set(gca,'XTick',TTick);
	set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
	set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
	xlim([DateNumBeg DateNumEnd]);

%% Calculating Y axis and Color axis limits
	miny=min(yAxis); maxy=max(yAxis);
	ylim([miny maxy]);

	minz=min(min(zValue)); maxz=max(max(zValue));
	caxis([minz maxz]);
	colorbar;

%% Defining necessary labels according to different modes
	switch mode
		case 1
			set(get(gca,'YLabel'),'String','Altitude [km]');
			title(['Electron Density log_1_0 [m^-^3] at ',datestr(time(1),'DD-MM-YY')]);
		case 2
			set(get(gca,'YLabel'),'String','Altitude [km]');
			title(['Production Rates log_1_0 [m^-^3 s^-^1] at ',datestr(time(1),'DD-MM-YY')]);
		case 3
			set(get(gca,'YLabel'),'String','Energy [eV]');
			title(['Differential Energy Flux log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1] at ',datestr(time(1),'DD-MM-YY')]);
			set(gca,'YScale','log');
		otherwise
		    warning('No mode selected')
	end



end