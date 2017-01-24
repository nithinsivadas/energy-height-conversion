function p=plot_1D_time_slice(time, yAxis, zValue, thisTime, mode)

%% plot2D.m Plots the electron energy flux using an input data structure

%  time 	: 1xN - A time array [matlab units]
%  yAxis 	: Mx1 - Values of the parameter on vertical axis (y-axis)
%                 Example: Energy bin [eV]
%  zValue 	: MxN - A matrix specifying zValue along time and y-axis
%                 Example: Energy flux [eV/ cm^2 sec sr eV]
%  timeTime : Value of time at which to plot [String]
%  mode     : Uses a set of pre-determined plot options
%       '-1': Assumes loglog(yAxis,zValueThisTime) with no labels [Default]
%        '0': Assumes semilogx(zValueThisTime,yAxis) with no labels 
%		 '1': Assumes one is plotting Electron Density
%        '2': Assumes one is plotting Production rates
%        '3': Assumes one is plotting Differential Energy Flux
%----------------------------------------------------------------------------
% Modified: 26th Sep 2016 : added plot handle 'p'
% Created : 22nd Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

%% Default Settings
	if nargin>5
	    error('plot_eflux:TooManyInputs','required at most 7 inputs');
    end;
    if nargin<5
        mode=-1; % Default Mode
    end;
    if nargin<4
	    error('plot_eflux:TooFewInputs','required at least 4 inputs');
	end;
    
    timeNo = find_time(time,thisTime);
    zValueThisTime = zValue(:,timeNo);
	grid on;
    
    switch mode
    	case 3
    		p=plot(yAxis,zValueThisTime);
	        set(gca,'XScale','log');
	        set(gca,'YScale','log');
	    	
	    	set(get(gca,'XLabel'),'String','Energy [eV]');
	    	set(get(gca,'YLabel'),'String','Differential Energy Flux [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
			legend(datestr(time(timeNo)));	

		case 2
	    	p=plot(zValueThisTime,yAxis);
	    	set(gca,'XScale','log');
			set(get(gca,'XLabel'),'String','Production Rates [m^-^3 s^-^1]');
	    	set(get(gca,'YLabel'),'String','Altitude [km]');			

	    case 1
	       	p=plot(zValueThisTime,yAxis);
    		set(gca,'XScale','log');
    		set(get(gca,'XLabel'),'String','Electron Density [m^-^3]');
    		set(get(gca,'YLabel'),'String','Altitude [km]');

    	case 0
	    	p=plot(zValueThisTime,yAxis);
	    	set(gca,'XScale','log');

    	otherwise
    		p=plot(yAxis,zValueThisTime);
    		set(gca,'XScale','log');
	        set(gca,'YScale','log');
    end

%     set(gca, 'Layer','top')


end

