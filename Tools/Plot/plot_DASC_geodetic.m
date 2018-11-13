function [ax1,h] = plot_DASC_geodetic( dataNew, time, lat, lon,...
    imageSize, latLim, lonLim, deltaLat, deltaLon, label )
%% plot_DASC_geodetic Plot DASC optical image in geodetic coordinates on a map
%--------------------------------------------------------------------------
% Input
%------
% dataNew - 2-D optical data from DASC in geodetic coordinates [nCoordinates]
%           produced from DASC_aer_to_geodetic.m
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
if nargin<8
    deltaLat = 1;
end

if nargin<9
    deltaLon = 5;
end

if nargin<10
    label = 'DASC';
elseif strcmp(label,'pokerFlat')
    label = 'PKFT DASC';
end
nanFlags = isnan(lat) | isnan(lon) | isnan(dataNew);
lat(nanFlags) = [];
lon(nanFlags) = [];
dataNew(nanFlags) = [];
% lat, lon, dataNew have to be column vectors!
F = scatteredInterpolant(lat',lon',dataNew','nearest','none'); %Modified on 3rd Oct 2018 - Transpose
latq = linspace(min(lat),max(lat),imageSize);
lonq = linspace(min(lon),max(lon),imageSize);
[LAT, LON] = ndgrid(latq,lonq);
Vq = F(LAT,LON);
ax1=axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);
axis off
load coastlines
plotm(coastlat,coastlon)
hold on;
h=pcolorm(LAT,LON,(Vq)); 
set(h,'EdgeColor','none');
if nargin>9
    textm(latLim(2)+(latLim(2)-latLim(1))*0.05, lonLim(1) +(lonLim(2)-lonLim(1))*0.35, [char(upper(label)),': ',datestr(time,'HH:MM:SS'),' UT']);
end
% textm(latLim(2)-0.19, lonLim(1)+0.1, ['DASC: ',datestr(time,'HH:MM:SS'),' UT']);
end

