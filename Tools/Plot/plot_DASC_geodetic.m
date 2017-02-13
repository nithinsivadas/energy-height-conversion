function [ax1,h] = plot_DASC_geodetic( dataNew, time, lat, lon, imageSize, latLim, lonLim )
%% plot_DASC_geodetic Plot DASC optical image in geodetic coordinates on a map
%--------------------------------------------------------------------------
% Input
%------
% dataNew - 2-D optical data from DASC in geodetic coordinates [nCoordinates]
% time    - Matlab time of the particular DASC image [nCoordinates]
% lat     - latitude coordinates [nCoordinates]
% lon     - longitude coordinates [nCoordinates]
% imageSize - the total pixel size of the output image e.g. 1024 
% latLim    - latitude limits e.g. [61 65]
% lonLim    - longitude limits e.g. [141.5 144.5]
%--------------------------------------------------------------------------
% Output
%-------
% ax1 - map axes
% h   - pcolor handle of color plot of the optical data
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------

F = scatteredInterpolant(lat,lon,dataNew,'nearest','none');
latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);
[LAT, LON] = ndgrid(latq,lonq);
Vq = F(LAT,LON);

ax1=axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',1,'MLineLocation',5);
axis off
load coastlines
plotm(coastlat,coastlon)
hold on;
h=pcolorm(LAT,LON,(Vq)); 
set(h,'EdgeColor','none');
textm(latLim(2)-0.2, lonLim(1)+0.1, ['DASC: ',datestr(time,'HH:MM:SS'),' UT']);

end

