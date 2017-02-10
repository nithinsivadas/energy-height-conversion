function [h2] = plot_2D_energy_slice_aer( data, aercoords, energyBin, nBeams, timeNo, altitude, energy, setMapOn )
%% plot_2D_energy_slice_aer.m Plot 2D energy flux map at a particular time, energy and projected at an altitude
%--------------------------------------------------------------------------
% Input:
%-------
% data          - arranged beam-wise
%    -energyFlux- contains inverted energy flux varying with energy bin for each beam
%    -time      - time array
% aercoords     - azimuth, elevation, range coordinates arranged non-beam wise
% energyBin     - vector containing energy bin values of the data
% nBeams        - number of beams
% timeNo        - index of time value in time array 
% altitude      - altitude value where image ought to be projected
% energy        - energy array 
% setMapOn      - true: plots map axes
%--------------------------------------------------------------------------
% Output:
%-------
% h2            - pcolor handle
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------

if nargin < 8
    setMapOn = true;
end


%% Generating lat, lon, h, energy coordinates

az = aercoords(:,1);
el = aercoords(:,2);
zUp = aercoords(:,4);
altitudeGrid = zUp(1:nBeams:end);
altitudeNo=find_altitude(altitudeGrid(:,1),altitude);

az1  = az(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));
el1 = el(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));
zUp1    = zUp(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));
x1 = (90-el1).*sind(az1);
y1 = (90-el1).*cosd(az1); 

az_new = repmat(az1,1,length(energyBin));
el_new = repmat(el1,1,length(energyBin));
zUp    = repmat(zUp1,1,length(energyBin));
zEnergyBin = repmat(energyBin,1,length(az_new));

az_new = reshape(az_new',[],1);
el_new = reshape(el_new',[],1);
zUp = reshape(zUp',[],1);
zEnergyBin = reshape(zEnergyBin,[],1);

%% Generating data slice
nEnergyBin = length(energyBin);
for iBeam=1:1:nBeams
    diffEnergyFlux(1+(iBeam-1)*nEnergyBin:1:(iBeam)*nEnergyBin,1) =...
        data(iBeam).energyFlux(:,timeNo); % Choosing energy and time slice
end;
x = (90-el_new).*sind(az_new(:));
y = (90-el_new).*cosd(az_new(:)); 

F = scatteredInterpolant(x, y, zEnergyBin, diffEnergyFlux,'nearest','none');
imageSize = 512;
xLim = [min(x(:)) max(x(:))];
yLim = [min(y(:)) max(y(:))];

xq = linspace(xLim(1),xLim(2),imageSize);
yq = linspace(yLim(1),yLim(2),imageSize);
Vq = F({xq,yq,energy*1000});

if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[floor(xLim(1)) ceil(xLim(2))],...
        'MapLonLimit',[floor(yLim(1)) ceil(yLim(2))],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',1,'MLineLocation',1);
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
end;

h2=pcolor(xq,yq,log10(Vq)'); 
set(h2,'EdgeColor','none');
% hold on;
scatter(x1,y1,'.','w');
% caxis([8 12]);
hold on;
text(xLim(2), yLim(2)+1, ['PFISR: ', num2str(energy),' keV'],'color','r');
text(xLim(2)-2, yLim(2)+1,...
    [datestr(data(1).time(timeNo),'HH:MM:SS'),' UT'],'color', 'r');
end

