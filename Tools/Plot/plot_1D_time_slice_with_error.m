function p=plot_1D_time_slice_with_error(time, yAxis, zValue, zError, thisTime, mode)

%% plot_1D_time_slice.m Plots the 1-D energy flux using an input data structure
%--------------------------------------------------------------------------
% Input
%------
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
%--------------------------------------------------------------------------
% Output
%-------
% p         : plot handle
%--------------------------------------------------------------------------
% Modified: 26th Sep 2016 : added plot handle 'p'
% Created : 22nd Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%--------------------------------------------------------------------------
%%

%% Default Settings
	if nargin>6
	    error('plot_eflux:TooManyInputs','required at most 6 inputs');
    end;
    if nargin<6
        mode=-1; % Default Mode
    end;
    if nargin<5
	    error('plot_eflux:TooFewInputs','required at least 5 inputs');
	end;
    
    timeNo = find_time(time,thisTime);
    zValueThisTime = zValue(:,timeNo);
    zErrorThisTime = zError(:,timeNo);
	grid on;
    
    switch mode
    	case 3
    		p=errorbar(yAxis,zValueThisTime,zErrorThisTime,'vertical');
	        set(gca,'XScale','log');
	        set(gca,'YScale','log');
	    	
	    	set(get(gca,'XLabel'),'String','Energy [eV]');
	    	set(get(gca,'YLabel'),'String','Differential Energy Flux [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
			legend(datestr(time(timeNo)));	

		case 2
	    	p=errorbar(zValueThisTime,yAxis,zErrorThisTime,'horizontal');
	    	set(gca,'XScale','log');
			set(get(gca,'XLabel'),'String','Production Rates [m^-^3 s^-^1]');
	    	set(get(gca,'YLabel'),'String','Altitude [km]');			

	    case 1
	       	p=errorbar(zValueThisTime,yAxis,zErrorThisTime,'horizontal');
    		set(gca,'XScale','log');
    		set(get(gca,'XLabel'),'String','Electron Density [m^-^3]');
    		set(get(gca,'YLabel'),'String','Altitude [km]');

    	case 0
	    	p=errorbar(zValueThisTime,yAxis,zErrorThisTime,'horizontal');
	    	set(gca,'XScale','log');

    	otherwise
    		p=errorbar(yAxis,zValueThisTime,zErrorThisTime,'vertical');
    		set(gca,'XScale','log');
	        set(gca,'YScale','log');
    end

%     set(gca, 'Layer','top')


end

