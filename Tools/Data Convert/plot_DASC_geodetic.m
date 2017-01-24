function [ax1,h] = plot_DASC_geodetic( dataNew, time, lat, lon, imageSize, latLim, lonLim )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

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

