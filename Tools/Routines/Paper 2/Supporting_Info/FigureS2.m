%% Supporting Info Figure S2
clear all;
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
%%
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
timeMinStr = '26 Mar 2008 08:30';
timeMaxStr = '26 Mar 2008 12:30';
%%
pfisr.eflux = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux'),[3 2 1]);
pfisr.time = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/time'),[2 1]);
pfisr.energyBin = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyBin'),[2 1]);
pfisr.altitude = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/alt'),[2 1]);
pfisr.Ne = permute(h5read(h5FileStr,'/inputData/Ne'),[3 2 1]);


%%
timeMinIndx = find_time(padData.thd.time,timeMinStr);
timeMaxIndx = find_time(padData.thd.time,timeMaxStr);
tRange = timeMinIndx:1:timeMaxIndx;
time = padData.thd.time(tRange);
eflux = padData.thd.lcDiffEfluxLi(tRange,:)*10^4;
eflux(eflux<=0) = 0.1;
energyBin = padData.thd.energyBin;
p=create_panels(figure, 'totalPanelNo',3);
q = p(1);

q(1).select();
plot_2D_time_series(time, energyBin, log10(eflux)', 0.5, 0);
ylim([10^3,10^6]);
set(gca,'YScale','log','YTick',[1,10,30,100,300].*10^3,'YTickLabel',{'1','10','30','100','300'});
ylabel({'Themis-D loss-cone','energy flux','[keV]'});

colormap(inferno);
caxis([8, 12]);
c=colorbar_thin();
set(c,'YTick',[8,9,10,11,12],'YTickLabel',{'10^8','10^9','10^1^0','10^1^1','10^1^2'});
ylabel(c,'[eV/m^2 sr s eV]');


timeMinIndx = find_time(pfisr.time,timeMinStr);
timeMaxIndx = find_time(pfisr.time,timeMaxStr);
tRange = timeMinIndx:1:timeMaxIndx;
time = pfisr.time(tRange);
eflux = squeeze(pfisr.eflux(tRange,13,:));
eflux(eflux<=0) = 0.1;
eflux = movmean(eflux,5,1);
energyBin = pfisr.energyBin;
q(2).select();
plot_2D_time_series(time, energyBin, log10(eflux)', 0.5, 0);
ylim([10^3,10^6]);
set(gca,'YScale','log');
ylabel({'PFISR \phi(E)','energy flux','[keV]'});
colormap(inferno);
caxis([8, 12]);
set(gca,'YTick',[1,10,30,100,300].*10^3,'YTickLabel',{'1','10','30','100','300'});
c=colorbar_thin();
set(c,'YTick',[8,9,10,11,12],'YTickLabel',{'10^8','10^9','10^1^0','10^1^1','10^1^2'});
ylabel(c,'[eV/m^2 sr s eV]');

timeMinIndx = find_time(pfisr.time,timeMinStr);
timeMaxIndx = find_time(pfisr.time,timeMaxStr);
tRange = timeMinIndx:1:timeMaxIndx;
time = pfisr.time(tRange);
Ne = squeeze(pfisr.Ne(tRange,13,:));
Ne(Ne<=0) = 0.1;
% Ne = movmean(Ne,5,1);
alt = pfisr.altitude;
q(3).select();
plot_2D_time_series(time, alt, log10(Ne)', 0.5, 0);
ylim([60,200]);
set(gca,'YScale','linear');
ylabel({'PFISR N_e','[km]'});
colormap(gca,magma);
caxis([8, 12]);
set(gca,'YTick',[60,80,100,120,200]);
c=colorbar_thin();
set(c,'YTick',[8,9,10,11,12],'YTickLabel',{'10^8','10^9','10^1^0','10^1^1','10^1^2'});
label_time_axis(time,true,0.5,timeMinStr,timeMaxStr);
ylabel(c,'[m^-^3]');
