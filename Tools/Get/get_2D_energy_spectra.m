function [dataInv,coords] = get_2D_energy_spectra(amisrData,energyBin,...
    timeMinStr,timeMaxStr,altLim,setTypeofCoords)
%get_2D_energy_spectra Converts electron density to energy spectra and
%prepares a struct variable which can be used to plot 2D energy flux slices
%   Input
%   setTypeofCoords = 'magnetic' or 'original'
%   energyBin [nE x 1] -- important

if nargin < 6
    setTypeofCoords = 'magnetic';
end
if nargin < 5
    altLim = [70 200];
end
if nargin < 4
    timeMaxStr = datestr(max(amisrData.time(1,:)));
end
if nargin < 3
    timeMinStr = datestr(min(amisrData.time(1,:)));
end
if nargin < 2
    energyBin = logspace(3,6,30)';
end

timeMinIndx = find_time(amisrData.time(1,:),timeMinStr);
timeMaxIndx = find_time(amisrData.time(1,:),timeMaxStr);
% PFISR Location
lat = amisrData.site.latitude;
lon = amisrData.site.longitude;
nBeams = amisrData.nBeams;
timeRange = timeMinIndx:1:timeMaxIndx;
time = amisrData.time(1,timeRange);
hWait = waitbar(0);
if strcmp(setTypeofCoords,'magnetic')
    for iBeam = 1:1:nBeams
        custom_waitbar(hWait,iBeam,nBeams,'Calculating 2D energy fluxes');
        minAltNo = find_altitude(amisrData.altitude(:,iBeam),altLim(1));
        maxAltNo = find_altitude(amisrData.altitude(:,iBeam),altLim(2));
        coordRange = minAltNo:1:maxAltNo;
        Ne = squeeze(amisrData.magElectronDensity(coordRange,iBeam,timeRange));
        alt = amisrData.altitude(coordRange,iBeam);
        q = get_production_rate(Ne,alt,time,2); 
        A = get_energy_dep_matrix(alt,energyBin,lat,lon,time(round(length(time)/2)));% Approximation - should iterate time
        dataInv(iBeam) = get_inverted_flux(q,time,alt,energyBin,A);
        coords(:,iBeam,:) = [amisrData.magCartCoords.xEast(coordRange,iBeam),amisrData.magCartCoords.yNorth(coordRange,iBeam),amisrData.magCartCoords.zUp(coordRange,iBeam)];
    end
elseif strcmp(setTypeofCoords,'original')
    for iBeam = 1:1:nBeams
        custom_waitbar(hWait,iBeam,nBeams,'Calculating 2D energy fluxes');
        minAltNo = find_altitude(amisrData.altitude(:,iBeam),altLim(1));
        maxAltNo = find_altitude(amisrData.altitude(:,iBeam),altLim(2));
        coordRange = minAltNo:1:maxAltNo;
        Ne = squeeze(amisrData.electronDensity(coordRange,iBeam,timeRange));
        alt = amisrData.altitude(coordRange,iBeam);
        q = get_production_rate(Ne,alt,time,2); 
        A = get_energy_dep_matrix(alt,energyBin,lat,lon,time(round(length(time)/2)));% Approximation - should iterate time
        dataInv(iBeam) = get_inverted_flux(q,time,alt,energyBin,A);
        coords(:,iBeam,:) = [amisrData.cartCoords.xEast(coordRange,iBeam),amisrData.cartCoords.yNorth(coordRange,iBeam),amisrData.cartCoords.zUp(coordRange,iBeam)];
    end
else
    error('setTypeofCoords unknown');
end
delete(hWait);
end

