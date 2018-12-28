%% Plot PFISR Energy Spectra for different dates
h52008 = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';
h52010 = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20101018.001_bc_15sec-energyFlux_v85km.h5';
%%
iBeam = 13;
pfisr08.ne=h5read(h52008,'/inputData/Ne');
pfisr08.energyFlux = h5read(h52008,'/energyFluxFromMaxEnt/energyFlux');
pfisr08.energyBin = h5read(h52008,'/energyFluxFromMaxEnt/energyBin');
pfisr08.coords = h5read(h52008,'/magneticFieldAlignedCoordinates/magGeodeticLatLonAlt');
pfisr08.time = h5read(h52008,'/energyFluxFromMaxEnt/time');

pfisr10.ne=h5read(h52010,'/inputData/Ne');
pfisr10.energyFlux = h5read(h52010,'/energyFluxFromMaxEnt/energyFlux');
pfisr10.energyBin = h5read(h52010,'/energyFluxFromMaxEnt/energyBin');
pfisr10.coords = h5read(h52010,'/magneticFieldAlignedCoordinates/magGeodeticLatLonAlt');
pfisr10.time = h5read(h52010,'/energyFluxFromMaxEnt/time');

%%
timeMinStr = '18 Oct 2010 07:00';
timeMaxStr = '18 Oct 2010 10:00';
h=figure;
p=create_panels(h,'totalPanelNo',4,'demargin',15);
p(1,1).select();
yAxis = squeeze(pfisr10.coords(:,iBeam,3));
zValue = squeeze(pfisr10.ne(:,iBeam,:)); 
zValue(zValue<0) = 10^4;
zValue = log10(zValue);
colormap(inferno);
plot_2D_time_series(pfisr10.time,yAxis,zValue,0.5,1,timeMinStr,timeMaxStr);
caxis([6,12]);
colorbar('Location','eastoutside');

p(1,2).select();
yAxis = squeeze(pfisr10.energyBin)/1000;
zValue = squeeze(pfisr10.energyFlux(:,iBeam,:)); 
zValue(zValue<0) = 10^4;
zValue = log10(zValue);
colormap(inferno);
plot_2D_time_series(pfisr10.time,yAxis,zValue,0.5,3,timeMinStr,timeMaxStr);
caxis([6,12]);
colorbar('Location','eastoutside');

timeMinStr = '26 Mar 2008 09:00';
timeMaxStr = '26 Mar 2008 12:00';
p(1,3).select();
yAxis = squeeze(pfisr08.coords(:,iBeam,3));
zValue = squeeze(pfisr08.ne(:,iBeam,:)); 
zValue(zValue<0) = 10^4;
zValue = log10(zValue);
colormap(inferno);
plot_2D_time_series(pfisr08.time,yAxis,zValue,0.5,1,timeMinStr,timeMaxStr);
caxis([6,12]);
colorbar('Location','eastoutside');

p(1,4).select();
yAxis = squeeze(pfisr08.energyBin)/1000;
zValue = squeeze(pfisr08.energyFlux(:,iBeam,:)); 
zValue(zValue<0) = 10^4;
zValue = log10(zValue);
colormap(inferno);
plot_2D_time_series(pfisr08.time,yAxis,zValue,0.5,3,timeMinStr,timeMaxStr);
caxis([6,10]);
colorbar('Location','eastoutside');

figure;
plot_1D_time_slice(
