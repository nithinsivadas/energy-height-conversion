function h=plot_energyflux_keogram( data, aercoords, energyBin, nBeams, timeMinStr, timeMaxStr, altitude, energy )
%% plot_energyflux_keogram.m Plot 2-D keogram of enegry flux at particular energy
%--------------------------------------------------------------------------
% Input
%------
% data      - 2-D optical data from DASC in geodetic coordinates [nCoordinates]
% aercoords - Azimuth, elevation, and range coordinates [nCoordinates x 3]
% energyBin - Array with values of energy bin
% nBeams    - Total number of beams
% timeMinStr - Lower limit of time in string
% timeMaxStr - Upper limit of time in string
% altitude   - The altitude at which the keogram is projected 
% energy     - The energy to be plotted
%--------------------------------------------------------------------------
% Output
%------
% h          - pcolor plot handle of plot_2D_time_series
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------

az = aercoords(:,1);
el = aercoords(:,2);
zUp = aercoords(:,4);
altitudeGrid = zUp(1:nBeams:end);
altitudeNo=find_altitude(altitudeGrid(:,1),altitude);

az1  = az(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));
el1 = el(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));
zUp1    = zUp(1+nBeams*(altitudeNo-1):1:nBeams*(altitudeNo));

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
timeMinNo = find_time(data(1).time,datenum(timeMinStr));
timeMaxNo = find_time(data(1).time,datenum(timeMaxStr));
thisTime=1;
for itime = timeMinNo:1:timeMaxNo
    for iBeam=1:1:nBeams
        diffEnergyFlux(1+(iBeam-1)*nEnergyBin:1:(iBeam)*nEnergyBin,1) =...
            data(iBeam).energyFlux(:,itime); % Choosing energy and time slice
    end;
    x = (90-el_new).*sind(az_new(:));
    y = (90-el_new).*cosd(az_new(:)); 

    F = scatteredInterpolant(x, y, zEnergyBin, diffEnergyFlux,'nearest','none');

    xq = zeros(1,256);
    yq = linspace(-90,90, 256);
    eq = energy*1000*ones(1,256);
    Vq(:,thisTime) = F(xq, yq, eq);
    Vtime(thisTime) = data(1).time(itime);
    thisTime=thisTime+1;
end;
    
h=plot_2D_time_series(Vtime, yq, log10(Vq), 0.5, 0, timeMinStr, timeMaxStr);
set(gca,'YScale','linear');
ylabel('Elevation [deg]');


end

