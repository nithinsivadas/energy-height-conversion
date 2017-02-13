function [h2] = plot_2D_energy_slice_geodetic( data, magcoords, energyBin, nBeams, timeNo, altitude, energy, setMapOn )
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
% magcoords     : arranged non-beam-wise [nCoordinates x 3]
% energyBin     : Energy bin values [nE x 1]
% nBeams        : Total number of beams
% timeNo        : Time number of the energy slice to be plotted
% altitude      : Altitude of projection of the energy slice
% energy        : Energy in keV of the differential energy flux to be plotted
% setMapOn      : True => Map axis on
%--------------------------------------------------------------------------
%Output
%------
% h2 - plot handle
%
%% 
%----------------------------------------------------------------------------
% Modified: 24th Jan 2017 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%
if nargin < 8
    setMapOn = true;
end
%% Generating lat, lon, h, energy coordinates

[lat,lon,zMagUp ] = magcoords_ned2geodetic( magcoords);

xMagEast  = magcoords(:,1);
yMagNorth = magcoords(:,2);
zMagUp    = magcoords(:,3);

altitudeGrid = zMagUp(1:1:nBeams);
altitudeNo=find_altitude(altitudeGrid(:,1),altitude);

lat1  = lat(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));
lon1 = lon(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));
zUp1    = zMagUp(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));

latitude = repmat(lat1,1,length(energyBin));
longitude = repmat(lon1,1,length(energyBin));
zUp    = repmat(zUp1,1,length(energyBin));
zEnergyBin = repmat(energyBin,1,length(latitude));

latitude = reshape(latitude',[],1);
longitude = reshape(longitude',[],1);
zUp = reshape(zUp',[],1);
zEnergyBin = reshape(zEnergyBin,[],1);

%% Generating data slice
nEnergyBin = length(energyBin);
for iBeam=1:1:nBeams
    diffEnergyFlux(1+(iBeam-1)*nEnergyBin:1:(iBeam)*nEnergyBin,1) =...
        data(iBeam).energyFlux(:,timeNo); % Choosing energy and time slice
end;

F = scatteredInterpolant(latitude, longitude, zEnergyBin, diffEnergyFlux,'nearest','none');
imageSize = 512;
latLim = [min(latitude(:)) max(latitude(:))];
lonLim = [min(longitude(:)) max(longitude(:))];

latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);
Vq = F({latq,lonq,energy*1000});

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
hold on;
textm(latLim(2), lonLim(2)+0.1, ['PFISR: ', num2str(energy),' keV'],'color','r');
textm(latLim(2)-0.2, lonLim(2)+0.1,...
    [datestr(data(1).time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');
end

