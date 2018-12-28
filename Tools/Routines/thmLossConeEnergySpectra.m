%% Plot themis loss-cone energy spectra
clear all;
%% Import data
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS89\lossConeFluxTHMADE_20080326_TS89_with_footpoints.mat');

%% Calculating L-shells
omniH5='G:\My Drive\Research\Projects\Data\omni.h5';
[maginput,magTime] = generate_maginput(omniH5,'26-Mar-2008 08:00','26-Mar-2008 13:00');

%% THA
kext = find_irbem_magFieldModelNo('TS89');
maginput = filter_irbem_maginput(kext,maginput);
time = padData.tha.time;
xGSE = padData.tha.XYZ_GSE;
tic
for iTime = 1:1:length(time)
    thisTime = time(iTime);
    thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
    padData.tha.Lm(iTime,1)=onera_desp_lib_make_lstar(kext,[0,0,0,0,0],3,thisTime,xGSE(iTime,1),xGSE(iTime,2),xGSE(iTime,3),thisMaginput);
end
toc

%% THD
kext = find_irbem_magFieldModelNo('TS89');
maginput = filter_irbem_maginput(kext,maginput);
time = padData.thd.time;
xGSE = padData.thd.XYZ_GSE;
tic
for iTime = 1:1:length(time)
    thisTime = time(iTime);
    thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
    padData.thd.Lm(iTime,1)=onera_desp_lib_make_lstar(kext,[0,0,0,0,0],3,thisTime,xGSE(iTime,1),xGSE(iTime,2),xGSE(iTime,3),thisMaginput);
end
toc
%% THE
kext = find_irbem_magFieldModelNo('TS89');
maginput = filter_irbem_maginput(kext,maginput);
time = padData.the.time;
xGSE = padData.the.XYZ_GSE;
tic
for iTime = 1:1:length(time)
    thisTime = time(iTime);
    thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
    padData.the.Lm(iTime,1)=onera_desp_lib_make_lstar(kext,[0,0,0,0,0],3,thisTime,xGSE(iTime,1),xGSE(iTime,2),xGSE(iTime,3),thisMaginput);
end
padData.the.Lm(end)=interp1(1:1:length(padData.the.time)-1,padData.the.Lm(1:end-1),length(padData.the.time),'linear','extrap');
toc
%% 
cLim = [6 12];
timeMinStr = '26-Mar-2008 08:00';
timeMaxStr = '26-Mar-2008 13:00';

h = figure;
p = create_panels(h,'totalPanelNo',3,'demargin',15);

colormap(inferno);
time = padData.tha.time;
mlat = padData.tha.mlatFoot;
mlt = padData.tha.mltFoot;
Lm = padData.tha.Lm;
yAxis = padData.tha.energyBin/1000; % Converting to keV
zValue = padData.tha.lcDiffEfluxLi' * 10^4; % Converting to eV/m2 sr s eV
zValue(zValue<0) = 10^3;
p(1,1).select();
title('THEMIS Loss cone flux : 26 Mar 2008');
plot_2D_time_series(time,yAxis,log10(zValue),0.5,0);
set(gca,'YTick',[1,10,100,300],'YTickLabel',{'1','10','100','300'},...
    'YLim',[1,1000],'YScale','log');
ylabel({'Themis-A','Loss-cone flux','[keV]'});
caxis(cLim);
cb1=colorbar('eastoutside','Ticks',[6,8,10,12],'TickLabels',{'10^6','10^8','10^1^0','10^1^2'});
cb1.Position(1) = 0.9;
cb1.Position(3) = 0.01;
cb1.Label.String = '[eV m^-^2 s^-^1 sr^-^1 eV^-^1]';
[TTick,TTickLim]=label_time_axis(time,false,0.5,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,time,mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,time,mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,time,Lm,'Lm',3);

time = padData.thd.time;
mlat = padData.thd.mlatFoot;
mlt = padData.thd.mltFoot;
Lm = padData.thd.Lm;
yAxis = padData.thd.energyBin/1000;
zValue = padData.thd.lcDiffEfluxLi' * 10^4;
zValue(zValue<0) = 10^3;
p(1,2).select();
plot_2D_time_series(time,yAxis,log10(zValue),0.5,0);
set(gca,'YTick',[1,10,100,300],'YTickLabel',{'1','10','100','300'},...
    'YLim',[1,1000],'YScale','log');
ylabel({'Themis-D','Loss-cone flux','[keV]'});
caxis(cLim);
cb2=colorbar('eastoutside','Ticks',[6,8,10,12],'TickLabels',{'10^6','10^8','10^1^0','10^1^2'});
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[eV m^-^2 s^-^1 sr^-^1 eV^-^1]';
[TTick,TTickLim] = label_time_axis(time,false,0.5,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,time,mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,time,mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,time,Lm,'Lm',3);
cumEnergyFlux = diff_to_cumu_flux(zValue(30:end,:),yAxis(30:end));
medianIndx = find_median_energy(cumEnergyFlux,1);
cumEnergyBin = yAxis(30:end);
plot(time,cumEnergyBin(medianIndx),'c');

cumEnergyFluxLow = diff_to_cumu_flux(zValue(1:30,:),yAxis(1:30));
medianIndx = find_median_energy(cumEnergyFluxLow,1);
cumEnergyBinLow = yAxis(1:30);
plot(time,cumEnergyBinLow(medianIndx),'m');

time = padData.the.time;
mlat = padData.the.mlatFoot;
mlt = padData.the.mltFoot;
Lm = padData.the.Lm;
yAxis = padData.the.energyBin/1000;
zValue = padData.the.lcDiffEfluxLi' * 10^4;
zValue(zValue<0) = 10^3;
p(1,3).select();
plot_2D_time_series(time,yAxis,log10(zValue),0.5,0);
set(gca,'YTick',[1,10,100,300],'YTickLabel',{'1','10','100','300'},...
    'YLim',[1,1000],'YScale','log');
ylabel({'Themis-E','Loss-cone flux','[keV]'});
caxis(cLim);
cb3=colorbar('eastoutside','Ticks',[6,8,10,12],'TickLabels',{'10^6','10^8','10^1^0','10^1^2'});
cb3.Position(1) = 0.9;
cb3.Position(3) = 0.01;
cb3.Label.String = '[eV m^-^2 s^-^1 sr^-^1 eV^-^1]';
[TTick,TTickLim] = label_time_axis(time,true,0.5,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,time,mlat,'MLAT',2);
add_horizontal_axes(TTick,TTickLim,time,mlt,'MLT',3);
add_horizontal_axes(TTick,TTickLim,time,Lm,'Lm',4);
cumEnergyFlux = diff_to_cumu_flux(zValue(30:end,:),yAxis(30:end));
medianIndx = find_median_energy(cumEnergyFlux,1);
cumEnergyBin = yAxis(30:end);
plot(time,cumEnergyBin(medianIndx),'c');

cumEnergyFluxLow = diff_to_cumu_flux(zValue(1:30,:),yAxis(1:30));
medianIndx = find_median_energy(cumEnergyFluxLow,1);
cumEnergyBinLow = yAxis(1:30);
plot(time,cumEnergyBinLow(medianIndx),'m');
