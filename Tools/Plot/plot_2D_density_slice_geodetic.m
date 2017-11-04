function [ figureHandle ] = plot_2D_density_slice_geodetic...
    ( data, geodeticCoords, timeNo, altitudeSelected, setMapOn, latLim, lonLim)
%plot_2D_density_slices_geodetic Plotting a 2D slice of PFISR densities
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
    setMapOn=true;
end

%% Initializing coordinates
lat = geodeticCoords(:,1);
lon = geodeticCoords(:,2);
alt = geodeticCoords(:,3);

electronDensity = data.electronDensity(:,timeNo);

%% Generating data slice
F = scatteredInterpolant(lat, lon, alt, electronDensity, 'nearest', 'none');
imageSize = 512;

% Loading default values for the lat and lon limit, if not specified
if nargin<6
latLim = [min(lat(:)) max(lat(:))];
end
if nargin<7
lonLim = [min(lon(:)) max(lon(:))];
end
%% Plotting map
latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);
Vq = F({latq,lonq,altitudeSelected});

if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[floor(latLim(1)) ceil(latLim(2))],...
        'MapLonLimit',[floor(lonLim(1)) ceil(lonLim(2))],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',1,'MLineLocation',1);
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
    h2=pcolorm(latq,lonq,log10(Vq)); 
    hold on;
    textm(latLim(2), lonLim(1)+0.1, ['PFISR: ', num2str(altitudeSelected),' km'],'color','r');
    textm(latLim(2)-0.2, lonLim(1)+0.1,...
     [datestr(data.time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');

else
    h2=pcolor(lonq,latq,log10(Vq)); 
    hold on;
    text(latLim(2), lonLim(2)+0.1, ['PFISR: ', num2str(altitudeSelected),' km'],'color','r');
    text(latLim(2)-0.2, lonLim(2)+0.1,...
     [datestr(data.time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');

end

set(h2,'EdgeColor','none');

%check plot_2D_energy_slice_geodetic.m
end

