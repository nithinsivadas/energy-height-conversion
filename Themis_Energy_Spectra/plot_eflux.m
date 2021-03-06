function [ ] = plot_eflux( data, time_tick, title_str )
%plot_eflux Plots the electron energy flux using an input data structure
%  data - 3 struct
%  data.time - Nx1 - A time array [matlab units]
%  data.ebin - 1xM - An array of energy values correspinding to each
%                    energy bin [eV]
%  data.E    - NxM - A matrix of energy flux [eV/ cm^2 sec sr eV]
%  time_tick - Value of time interval between tick marks [Hr]
%  title_str - Title string

% 12 May 2016
% Nithin Sivadas


if nargin>3
    error('plot_eflux:TooManyInputs','required at most 2 inputs');
end;

switch nargin
    case 1
        time_tick=0.5;
        title_str = 'Electron Energy Flux in log10 [eV/ cm^2 sec sr eV]';
        display('Setting default value of 0.5 hr spaced time tick points');
    case 2
        title_str = 'Electron Energy Flux in log10 [eV/ cm^2 sec sr eV]';
end

time = data.time;
ebin = data.ebin;
E    = data.E;

cmin=2; cmax=9;
minz=min(ebin); maxz=max(ebin);

[X,Y]=meshgrid(time,ebin);	
surf(X,Y,real(log10((E)')),'EdgeColor','none');

DateNumBeg=min(time);
DateNumEnd=max(time);

% Creating Time tick marks
dt=time_tick;
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
set(gca, 'Layer','top') 
set(gca,'XTick',TTick);
set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Energy [eV]');
title(title_str);
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax]);colorbar;


end

