%% Paper 2: Figure 3a-c
% Plot themis loss-cone energy spectra and L-shells
clear all;
%% Import data
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
% Load POES
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);

%% Calculating L-shells
omniH5='G:\My Drive\Research\Projects\Data\omni.h5';
[maginput,magTime] = generate_maginput(omniH5,'26-Mar-2008 08:00','26-Mar-2008 13:00');

%% Calculating Flux perpendicular to the field line
for iE = 1:1:size(padData.thd.fitAr.x0,2)
    padData.thd.trappedDiffEfluxAr(:,iE) = padData.thd.pitchAngleDistributionFunctionAr(squeeze(padData.thd.fitAr.x0(:,iE,:)),90);
    padData.the.trappedDiffEfluxAr(:,iE) = padData.the.pitchAngleDistributionFunctionAr(squeeze(padData.the.fitAr.x0(:,iE,:)),90);
    padData.thd.trappedDiffEfluxYi(:,iE) = padData.thd.pitchAngleDistributionFunctionYi(squeeze(padData.thd.fitYi.x0(:,iE,:)),90);
    padData.the.trappedDiffEfluxYi(:,iE) = padData.the.pitchAngleDistributionFunctionYi(squeeze(padData.the.fitYi.x0(:,iE,:)),90);
end

%%
smoothFactor = 10;
kext = find_irbem_magFieldModelNo('TS96');
maginput = filter_irbem_maginput(kext,maginput);

% POES trajectory
poesTimeMinStr = '26-Mar-2008 11:27';
poesTimeMaxStr = '26-Mar-2008 11:33';
poesTimeIndx = find_time(poes.time,poesTimeMinStr):1:find_time(poes.time,poesTimeMaxStr);
poestime = poes.time(poesTimeIndx);
for i = 1:1:length(poestime)
    poesGSM(i,:) = onera_desp_lib_coord_trans([850,poes.lat(poesTimeIndx(i)),poes.lon(poesTimeIndx(i))],...
        [0 2],poestime(i));
    thisMaginput = interp1(magTime,maginput,poestime(i));
    poesLm(i,1) = onera_desp_lib_make_lstar(kext,[0,0,0,0,0],2,poestime(i),...
        poesGSM(i,1),poesGSM(i,2),poesGSM(i,3),thisMaginput);
end

% THD
time = padData.thd.time;
xGSE = padData.thd.XYZ_GSE;
tic
for iTime = 1:1:length(time)
    thisTime = time(iTime);
    thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
    padData.thd.Lm(iTime,1)=onera_desp_lib_make_lstar(kext,[0,0,0,0,0],3,thisTime,xGSE(iTime,1),xGSE(iTime,2),xGSE(iTime,3),thisMaginput);
end
padData.thd.Lm = interp_nans(padData.thd.Lm);
padData.thd.LmSmooth = conv(padData.thd.Lm,ones(smoothFactor,1)/smoothFactor,'same');
toc

% THE

time = padData.the.time;
xGSE = padData.the.XYZ_GSE;
tic
for iTime = 1:1:length(time)
    thisTime = time(iTime);
    thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
    padData.the.Lm(iTime,1)=onera_desp_lib_make_lstar(kext,[0,0,0,0,0],3,thisTime,xGSE(iTime,1),xGSE(iTime,2),xGSE(iTime,3),thisMaginput);
end
padData.the.Lm(end)=interp1(1:1:length(padData.the.time)-1,padData.the.Lm(1:end-1),length(padData.the.time),'linear','extrap');
padData.the.Lm = interp_nans(padData.the.Lm);
padData.the.LmSmooth = conv(padData.the.Lm,ones(smoothFactor,1)/smoothFactor,'same');
toc


%%
cLim = [6 12];
timeMinStr = '26-Mar-2008 09:30';
timeMaxStr = '26-Mar-2008 12:30';
midEnergy = 30; % in KeV
h = figure;
p = create_panels(h,'totalPanelNo',3,'demargin',15);

% L-Shell variations
p(1,1).select();
scatter(padData.thd.time,padData.thd.Lm,5,[0.5,0.5,1],'filled');
hold on;
scatter(padData.the.time,padData.the.Lm,5,[0.5,1,1],'filled');
scatter(poestime,poesLm,10,[1,0.5,0.5],'filled');
p1=plot(padData.thd.time,padData.thd.LmSmooth,'b');
p2=plot(padData.the.time,padData.the.LmSmooth,'c');
p3=plot(poestime,poesLm,'r');
xlim([datenum(timeMinStr), datenum(timeMaxStr)]);
set(gca,'XTickLabel','');
ylabel({'TS96 L_m','[R_E]'});
[TTick,TTickLim] = label_time_axis(padData.the.time,true,0.5,timeMinStr,timeMaxStr);
ylim([8,15]);
legend([p1 p2 p3], {'Themis-D','Themis-E','NOAA-17'});

% THD
colormap(inferno);
time = padData.thd.time;
mlat = conv(padData.thd.mlatFoot,ones(smoothFactor,1)/smoothFactor,'same');
mlt = conv(padData.thd.mltFoot,ones(smoothFactor,1)/smoothFactor,'same');
Lm = padData.thd.LmSmooth;
yAxis = padData.thd.energyBin/1000;
zValue = padData.thd.lcDiffEfluxAr' * 10^4;
zValue(zValue<0) = 10^3;
% zValue1 = padData.thd.lcDiffEfluxAr'./padData.thd.trappedDiffEfluxAr';
midEBin = yAxis(find_altitude(yAxis,midEnergy));

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
cumEnergyFlux = diff_to_cumu_flux(zValue(midEBin:end,:),yAxis(midEBin:end));
medianIndx = find_median_energy(cumEnergyFlux,1);
cumEnergyBin = yAxis(midEBin:end);
plot(time,cumEnergyBin(medianIndx),'c');

cumEnergyFluxLow = diff_to_cumu_flux(zValue(1:midEBin,:),yAxis(1:midEBin));
medianIndx = find_median_energy(cumEnergyFluxLow,1);
cumEnergyBinLow = yAxis(1:midEBin);
plot(time,cumEnergyBinLow(medianIndx),'m');

% THE
colormap(inferno);
time = padData.the.time;
mlat = conv(padData.the.mlatFoot,ones(smoothFactor,1)/smoothFactor,'same');
mlt = conv(padData.the.mltFoot,ones(smoothFactor,1)/smoothFactor,'same');
Lm = padData.the.Lm;
yAxis = padData.the.energyBin/1000;
zValue = padData.the.lcDiffEfluxAr' * 10^4;
zValue(zValue<0) = 10^3;
% zValue1 = padData.the.lcDiffEfluxAr'./padData.the.trappedDiffEfluxAr';
midEBin = yAxis(find_altitude(yAxis,midEnergy));

p(1,3).select();
plot_2D_time_series(time,yAxis,log10(zValue),0.5,0);
% plot_2D_time_series(time,yAxis,zValue1,0.5,0);
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
cumEnergyFlux = diff_to_cumu_flux(zValue(midEBin:end,:),yAxis(midEBin:end));
medianIndx = find_median_energy(cumEnergyFlux,1);
cumEnergyBin = yAxis(midEBin:end);
plot(time,cumEnergyBin(medianIndx),'c');

cumEnergyFluxLow = diff_to_cumu_flux(zValue(1:midEBin,:),yAxis(1:midEBin));
medianIndx = find_median_energy(cumEnergyFluxLow,1);
cumEnergyBinLow = yAxis(1:midEBin);
plot(time,cumEnergyBinLow(medianIndx),'m');
