function [ figureHandle ] = plot_2D_density_slice_geodetic_aer...
    ( data, geodeticCoords, altitude, nBeams, timeNo, altitudeSelected, setMapOn)
% plot_2D_density_slice_geodetic_aer Plotting 2D slice in geodetic of density profiles
%   Detailed explanation goes here

if nargin <7
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
latLim = [min(lat(:)) max(lat(:))];
lonLim = [min(lon(:)) max(lon(:))];

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
end;

h2=pcolorm(latq,lonq,log10(Vq));
set(h2,'EdgeColor','none');
% hold on;
% textm(latLim(2), lonLim(2)+0.1, ['PFISR: ', num2str(energy),' keV'],'color','r');
% textm(latLim(2)-0.2, lonLim(2)+0.1,...
%     [datestr(data(1).time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');

%check plot_2D_energy_slice_geodetic.m
end
