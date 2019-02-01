function [ax2, h2] = plot_2D_energy_slice_geodetic_v2018( diffEnergyFlux,...
    latitude, longitude, zEnergyBin, timeNumPFISR, ...
    thisEnergy, latWidth, lonWidth, setMapOn,...
    setTimeLabelOn, imageSize)
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
if nargin < 11
    imageSize = 64; %Hard coded? based on cal file
end
if nargin < 10
    setTimeLabelOn = true;
end
if nargin < 9
    setMapOn = true;
end
if nargin < 8
    lonWidth = 4;
end
if nargin < 7
    latWidth = 2;
end
%% Generating lat, lon, h, energy coordinates
% diffEnergyFlux(energy,beams) - at a time instant
% lat, lon, energyBin, diffenergyflux - for all data points
%% Generating data slice
F = scatteredInterpolant(latitude(:), longitude(:), zEnergyBin(:), diffEnergyFlux(:),'nearest','none');

latLim = [min(latitude(:)) max(latitude(:))];
lonLim = [min(longitude(:)) max(longitude(:))];

latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);

Vq = F({latq,lonq,thisEnergy*1000});
Vq(Vq<=0)=nan;

if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1))-latWidth/2 (latLim(2))+latWidth/2],...
        'MapLonLimit',[(lonLim(1))-lonWidth/2 (lonLim(2))+lonWidth/2],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',1,'MLineLocation',1);
    
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
    h2=pcolorm(latq,lonq,log10(Vq)); 
else
    h2=pcolor(latq',lonq',log10(Vq)'); 
    
end


set(h2,'EdgeColor','none');

if setTimeLabelOn==true
    hold on;
    textm(latLim(2), lonLim(2)+lonWidth/10, ['PFISR: ', num2str(thisEnergy),' keV'],'color','r');
    textm(latLim(2)-latWidth/20, lonLim(2)+lonWidth/10,...
        [datestr(timeNumPFISR,'HH:MM:SS'),' UT'],'color', 'r');
    hold off;
end

end