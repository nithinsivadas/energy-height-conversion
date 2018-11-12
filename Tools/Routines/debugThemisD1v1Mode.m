%% Debug ThemisD1.v1 PFISR Mode
clear all;
%%
% inputH5Str = 'G:\My Drive\Research\Projects\Data\ThemisD1.v01\20180911.001_bc_1min-fitcal.h5';
inputH5Str = 'G:\My Drive\Research\Projects\Data\ThemisD1.v01\20180506.003_bc_1min-fitcal.h5';
 
% ouputH5Str = 'G:\My Drive\Research\Projects\Data\ThemisD1.v01\20180911.001_bc_1min-energyFlux.h5';

minAlt = 50;
maxAlt = 120;
projectionAltPFISR = 50; % km, altitude of origin of magnetic field aligned lines
nEnergyBins = 30;
minE = 5*10^3;
maxE = 10^6;

energyBin = logspace(minE,maxE,nEnergyBins);

%% Check Coordinates

pfisr = read_amisr(inputH5Str);
pfisrData = aer_to_field_aligned_coords(pfisr, projectionAltPFISR);

%% Check interpolation
minTimeStr = '11 Sep 2018 06:00:55';
maxTimeStr = '11 Sep 2018 06:10';
interpData = interpolate_to_field_aligned_coords(pfisrData,minTimeStr,maxTimeStr);
% There are some nan values on beam 5 consistently, suggesting this has to
% do something with the coordinates. 
%% Check Energy Inversion
[dataInv, magcoords, dataInputInv] = get_2D_energy_spectra(interpData,energyBin',...
   minTimeStr,maxTimeStr,[minAlt maxAlt],'magnetic',datenum('11-Sep-2008 06:00'));
%%
figure;
plot3(interpData.origCartCoords.xEast,interpData.origCartCoords.yNorth,...
    interpData.origCartCoords.zUp,'.-k');
hold on;
plot3(interpData.magCartCoords.xEast,interpData.magCartCoords.yNorth,...
    interpData.magCartCoords.zUp,'r');

%%
figure;
iBeam = 13;
minAlt = 30;
maxAlt = 120;
timeMinStr = '06 May 2018 15:00';
timeMaxStr = '06 May 2018 16:00';
altRangeIndx= find_altitude(pfisrData.altitude,minAlt):1:find_altitude(pfisrData.altitude,maxAlt);
figure;
zValue = (squeeze(pfisrData.electronDensity(altRangeIndx,iBeam,:)));
zValue(isnan(zValue))=10^6;
zValue(zValue<0)=10^6;
% zValue(~isreal(zValue)) = 1;
plot_2D_time_series(pfisrData.time(1,:),pfisrData.altitude(altRangeIndx,iBeam),log10(zValue),0.5,1,timeMinStr,timeMaxStr);
colormap('inferno');
colorbar;
caxis([10 12]);

%%
figure;

figure;
altRangeIndx= find_altitude(interpData.altitude,minAlt):1:find_altitude(interpData.altitude,maxAlt);
timeIndx = 10;
% beamIndx = interpData.magBeamNo;
beamIndx = 5;
plot(isnan(interpData.magElectronDensity(altRangeIndx,beamIndx,timeIndx)),...
    interpData.magGeodeticCoords.alt(altRangeIndx,beamIndx));
% set(gca,'XScale','log');
