function [h2] = plot_2D_energy_slice_geodetic_v1( data, amisrData, coords,...
    energyBin, nBeams, timeNo, altitude, energy, latWidth, lonWidth, setMapOn,...
    setTimeLabelOn)
%plot_2D_energy_slice_geodetic.m Plot 2D differential energy flux slices 
%from 4-D PFISR data sets on lat, long map
%--------------------------------------------------------------------------
%Input
%-----
% data          : arranged beam-wise
% -> flux       : differential number flux [nE x nTime]
% -> energyFlux : differential energy flux [nE x nTime]
% -> chi2       : Reduced chi-2 of the maximum entropy regression [1 x nTime]
% -> qInvert    : Production rate from inverted energy flux [nh x nTime]
% -> maxIter    : Maximum number of iterations [1 x nTime]
% -> energyBin  : Energy bin values [nE x 1]
% -> time       : Time array [nTime x 1]
% -> alt        : Altitude array [nh x 1]
% -> A          : Production rate vs. number flux matrix [nh x nE]
% -> qInput     : Production rate derived from electron density [nh x nTime]
% amisrData     : 
% -> site       :latitude, longitude and altitude (PFISR location)
% -> magBeamNo  :the beam number/ID that points along the mag. field line 
% magcoords     : arranged non-beam-wise [nh x 3 x nBeams]
% energyBin     : Energy bin values [nE x 1]
% nBeams        : Total number of beams
% timeNo        : Time number of the energy slice to be plotted
% altitude      : Altitude of projection of the energy slice
% energy        : Energy in keV of the differential energy flux to be plotted
% setMapOn      : True => Map axis on
% setTimeLabelOn: True => The time and energy values are printed on the
%                 plot
%--------------------------------------------------------------------------
%Output
%------
% h2 - plot handle
%
%% 
%----------------------------------------------------------------------------
% Modified: 2nd Feb 2018 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%
if nargin < 12
    setTimeLabelOn = true;
end
if nargin < 11
    setMapOn = true;
end
if nargin < 10
    lonWidth = 4;
end
if nargin < 9
    latWidth = 2;
end
%% Generating lat, lon, h, energy coordinates
lat0 = amisrData.site.latitude;
lon0 = amisrData.site.longitude;
h0 = amisrData.site.altitude/1000; %in Km
xMagEast  = coords(:,:,1);
yMagNorth = coords(:,:,2);
zMagUp    = coords(:,:,3);

[lat,lon] = ned2geodetic(yMagNorth,xMagEast,-zMagUp, lat0,lon0,h0,wgs84Ellipsoid('km'));

altitudeGrid = zMagUp(:,amisrData.magBeamNo);
altitudeNo=find_altitude(altitudeGrid,altitude);

lat1  = lat(altitudeNo,:);
lon1 = lon(altitudeNo,:);
zUp1    = zMagUp(altitudeNo,:);

latitude = repmat(lat1,length(energyBin),1);
longitude = repmat(lon1,length(energyBin),1);
zUp    = repmat(zUp1,length(energyBin),1);
zEnergyBin = repmat(energyBin,1,size(latitude,2));

%% Generating data slice
for iBeam=1:1:nBeams
    diffEnergyFlux(:,iBeam) =...
        data(iBeam).energyFlux(:,timeNo); % Choosing energy and time slice
end

F = scatteredInterpolant(latitude(:), longitude(:), zEnergyBin(:), diffEnergyFlux(:),'nearest','none');
imageSize = 512;
latLim = [min(latitude(:)) max(latitude(:))];
lonLim = [min(longitude(:)) max(longitude(:))];

latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);
Vq = F({latq,lonq,energy*1000});

if setMapOn==true
%     ax2=axesm('lambertstd','MapLatLimit',[floor(latLim(1)) ceil(latLim(2))],...
%         'MapLonLimit',[floor(lonLim(1)) ceil(lonLim(2))],...
%         'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
%         'PLineLocation',1,'MLineLocation',1);
    
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1))-latWidth/2 (latLim(2))+latWidth/2],...
        'MapLonLimit',[(lonLim(1))-lonWidth/2 (lonLim(2))+lonWidth/2],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',1,'MLineLocation',1);
    
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
end

h2=pcolorm(latq,lonq,log10(Vq)); 
set(h2,'EdgeColor','none');

% textm(latLim(2), lonLim(2)+0.1, ['PFISR: ', num2str(energy),' keV'],'color','r');
% textm(latLim(2)-0.2, lonLim(2)+0.1,...
%     [datestr(data(1).time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');
if setTimeLabelOn==true
    hold on;
    textm(latLim(2), lonLim(2)+lonWidth/10, ['PFISR: ', num2str(energy),' keV'],'color','r');
    textm(latLim(2)-latWidth/20, lonLim(2)+lonWidth/10,...
        [datestr(data(1).time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');
    hold off;
end
% plotm(lat1,lon1,'o');
end

