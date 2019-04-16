%% Event 2: 18 Oct 2010

%% Load PFISR Data
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20101018.001_bc_15sec-full_v110.h5';
data.Ne = h5read(h5FileStr,'/inputData/Ne'); % in [m^-3]
data.altitude = h5read(h5FileStr,'/energyFluxFromMaxEnt/alt');
data.time = h5read(h5FileStr,'/energyFluxFromMaxEnt/time');
%%
data.energyFlux = h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux');
data.energyBin = h5read(h5FileStr,'/energyFluxFromMaxEnt/energyBin');
%%
h= figure;
timeMinStr = '18 Oct 2010 07:30';
timeMaxStr = '18 Oct 2010 13:00';
% resize_figure(h,
iB =7;
electronDensity = data.Ne;
electronDensity(electronDensity<=0) = 1;
plot_2D_time_series(data.time,data.altitude,...
    log10(squeeze(electronDensity(:,iB,:))),1,0,timeMinStr, timeMaxStr);
colormap(inferno);
c=colorbar_thin();
ylabel(c,['log_1_0(N_e) Zenith Beam (No.',num2str(iB),')']);
ylabel({'PFISR N_e Densities','Altitude [km]'});
% label_time_axis(data.time(1,:),false,1,timeMinStr,timeMaxStr);
ylim([80 120]);
caxis([8 12]);
grid on;
title(datestr(floor(datenum(timeMinStr))));
label_time_axis(data.time,true,0.5,timeMinStr,timeMaxStr);

%%
figure;
energyFlux = data.energyFlux;
energyFlux(energyFlux<=0) = 1;
plot_2D_time_series(data.time,data.energyBin./1000,...
    log10(squeeze(energyFlux(:,iB,:))),1,0,timeMinStr, timeMaxStr);
colormap(inferno);
c=colorbar_thin();
ylabel(c,['log_1_0(\phi_E) Zenith Beam (No.',num2str(iB),')']);
ylabel({'PFISR Energy Spectra','Energy [keV]'});
% label_time_axis(data.time(1,:),false,1,timeMinStr,timeMaxStr);
ylim([1 1000]);
set(gca,'YScale','log');
caxis([8 10]);
grid on;
title(datestr(floor(datenum(timeMinStr))));
label_time_axis(data.time,true,0.5,timeMinStr,timeMaxStr);
