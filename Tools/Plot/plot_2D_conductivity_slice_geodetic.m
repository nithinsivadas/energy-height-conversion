function [ figureHandle ] = plot_2D_conductivity_slice_geodetic...
    ( data, geodeticCoords, timeNo, altitudeSelected, conductivityType, setMapOn, latLim, lonLim)
%plot_2D_conductivity_slice_geodetic Plotting a 2D slice of PFISR densities
% Input
% data : data.electronDensity --> Matrix of electron density point at each[position, time]
%      : the measurements can be organized in any order
%      : data.time --> Time index number
% geodeticCoords : [time x position] , position=(lat, lon, alt)
% timeNo : matlab datenum time of the map
% altitudeSelected: altitude slice of the map [km]
% setMapOn: 'true' value projects on a map
% latLim : the min and max value of the latitude to plot
% lonLim : the min and max value of the longitude to plot

%% Setting default values
if nargin <5
    conductivityType='H';
end

if nargin <6
    setMapOn=true;
end

%% Initializing coordinates
lat = geodeticCoords(:,1);
lon = geodeticCoords(:,2);
alt = geodeticCoords(:,3);

hallConductivity = data.sigma_H(:,timeNo);
pedersenConductivity = data.sigma_P(:,timeNo);

%% Generating data slice
if conductivityType=='H'
  F = scatteredInterpolant(lat, lon, alt, hallConductivity, 'nearest', 'none');
elseif conductivityType=='P'
  F = scatteredInterpolant(lat, lon, alt, pedersenConductivity, 'nearest', 'none');
end

imageSize = 512;

% Loading default values for the lat and lon limit, if not specified
if nargin<7
latLim = [min(lat(:)) max(lat(:))];
end
if nargin<8
lonLim = [min(lon(:)) max(lon(:))];
end
%% Plotting map
latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);
Vq = F({latq,lonq,altitudeSelected});

if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1)) (latLim(2))],...
        'MapLonLimit',[(lonLim(1)) (lonLim(2))],'FontSize',6,...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',0.5,'MLineLocation',1,'MLabelRound',-1,'PLabelRound',-1);
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
    h2=pcolorm(latq,lonq,log10(Vq));
    hold on;
%     titlem(latLim(2), lonLim(1)+0.1, [num2str(altitudeSelected),' km'],'color','r');
%     textm(latLim(2)-0.1, lonLim(1)+0.1,...
%      [datestr(data.time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');

else
    h2=pcolor(lonq,latq,log10(Vq));
    hold on;
    title(['PFISR: ', num2str(altitudeSelected),' km ',datestr(data.time(timeNo),'HH:MM:SS'),' UT']);
end

set(h2,'EdgeColor','none');

end
