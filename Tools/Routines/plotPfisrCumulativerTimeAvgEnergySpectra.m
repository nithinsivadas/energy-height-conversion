%%  Routine to Plot Cumulative Plots Energy Flux averaged across time
% We use SIC model, and constant alpha model to calculate q
% And the maximum entrop method to invert the spectra

clear all;
initialize_geodata;

%Global Variables
pfisrDataFileNameStr =...
    '/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/DataFile_2008_1.h5';
timeMin=datenum('26 Mar 2008 08:00');
timeMax=datenum('26 Mar 2008 13:00');
altMin = 60;
altMax = 150;

thmDataFileNameStr =...
    '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Espectra_thm.mat';
load (thmDataFileNameStr);
% data_thm -> E, ebin, time
thm.energyFlux = data_thm.E(:,19:1:end)'*10^4; %[eV m-2 sr-1 s-1 eV-1]
thm.energyBin  = data_thm.ebin(19:1:end)';
thm.time       = data_thm.time;
clear data_thm;
energyBin = thm.energyBin;

%%  Mode 1: From measurements averaged across the PFISR beams
[ePopDenAvg, altPop, timePop] = read_pfisr_variable(pfisrDataFileNameStr, 'popl', 0);

% Tidying up Electron Density
[NeAvg, time] = time_crop(ePopDenAvg, timePop, timeMin, timeMax);
[NeAvg, alt] = altitude_crop(NeAvg, altPop,altMin, altMax);
NeAvg = interp_nans(NeAvg);
NeAvg(NeAvg<0)=10^9;

% Generate Production Rate
[qAvg,qAvgTime, alpha] = get_production_rate(NeAvg, alt, time, 2);

% Invert production rates
[dataAvg] = get_inverted_flux(qAvg, qAvgTime, alt, energyBin);

%%  Mode 1.1: Error in PFISR Measurements
[dePopDenAvg, altPop, timePop] = read_pfisr_variable(pfisrDataFileNameStr, 'dpopl', 0);

% Tidying up Electron Density
[dNeAvg, time] = time_crop(dePopDenAvg, timePop, timeMin, timeMax);
[dNeAvg, alt] = altitude_crop(dNeAvg, altPop,altMin, altMax);
dNeAvg = interp_nans(dNeAvg);
dNeAvg(dNeAvg<0)=10^9;

% Generate Production Rate
[dqAvg,dqAvgTime, alpha] = get_production_rate(dNeAvg, alt, time, 2);

% Invert production rates
[dataAvgError] = get_inverted_flux(dqAvg, dqAvgTime, alt, energyBin);
%% Mode 2: From SIC Measurements

sicDataFileNameStr='/home/nithin/Documents/git-repos/ISR_School_2016/Antti  SIC Model/summaryresQinversion.mat';
load (sicDataFileNameStr);
[ionrate_raw, sic.alt] = altitude_crop(ionrate_raw, h,altMin, altMax);
[ionrate_smoothed, sic.alt] = altitude_crop(ionrate_smoothed, h,altMin, altMax);
sic.qRaw = ionrate_raw*10^6;
sic.qSmooth = ionrate_smoothed*10^6;
sic.Ne = Nedata*10^6;
sic.NeSIC = NeSIC*10^6;
sic.time = t;
clear h; clear ionrate_raw; clear ionrate_smoothed; clear Nedata; clear NeSIC; clear t; clear UT;

% Invert production rates
[dataAvgSIC] = get_inverted_flux(sic.qSmooth, sic.time, sic.alt, energyBin);

%% Calculating cumulative energy flux

dataAvg.cumEnergyFlux = ...
    diff_to_cumu_flux(dataAvg.energyFlux, dataAvg.energyBin);
dataAvgError.cumEnergyFlux = ...
    diff_to_cumu_flux(dataAvgError.energyFlux, dataAvgError.energyBin);
dataAvgSIC.cumEnergyFlux = ...
    diff_to_cumu_flux(dataAvgSIC.energyFlux, dataAvgSIC.energyBin);
thm.cumEnergyFlux = ...
    diff_to_cumu_flux(thm.energyFlux, thm.energyBin);

beforeOnset.timeMin = datenum('26-Mar-2008 11:00');
beforeOnset.timeMax = datenum('26-Mar-2008 11:44');

beforeOnset.dataAvg.cumuEnergyFlux = ...
    get_time_avg_time_series_data(dataAvg.cumEnergyFlux, dataAvg.time,...
    beforeOnset.timeMin, beforeOnset.timeMax);
beforeOnset.dataAvgError.cumuEnergyFlux = ...
    get_time_avg_time_series_data(dataAvgError.cumEnergyFlux,...
    dataAvgError.time,beforeOnset.timeMin, beforeOnset.timeMax);
beforeOnset.dataAvgSIC.cumuEnergyFlux = ...
    get_time_avg_time_series_data(dataAvgSIC.cumEnergyFlux,...
    dataAvgSIC.time,beforeOnset.timeMin, beforeOnset.timeMax);
beforeOnset.thm.cumuEnergyFlux = ...
    get_time_avg_time_series_data(thm.cumEnergyFlux,...
    thm.time,beforeOnset.timeMin, beforeOnset.timeMax);

afterOnset.timeMin = datenum('26-Mar-2008 11:44');
afterOnset.timeMax = datenum('26-Mar-2008 12:00');

afterOnset.dataAvg.cumuEnergyFlux = ...
    get_time_avg_time_series_data(dataAvg.cumEnergyFlux, dataAvg.time,...
    afterOnset.timeMin, afterOnset.timeMax);
afterOnset.dataAvgError.cumuEnergyFlux = ...
    get_time_avg_time_series_data(dataAvgError.cumEnergyFlux,...
    dataAvgError.time,afterOnset.timeMin, afterOnset.timeMax);
afterOnset.dataAvgSIC.cumuEnergyFlux = ...
    get_time_avg_time_series_data(dataAvgSIC.cumEnergyFlux,...
    dataAvgSIC.time,afterOnset.timeMin, afterOnset.timeMax);
afterOnset.thm.cumuEnergyFlux = ...
    get_time_avg_time_series_data(thm.cumEnergyFlux,...
    thm.time,afterOnset.timeMin, afterOnset.timeMax);

%% Creating Panels
hFig=figure(4);
resize_figure(hFig);

clf
p=panel(hFig);
% p.pack(1,2);
panelSize = 60; %in mm
p.pack({{panelSize} {panelSize}},{{80}});
p.marginleft=25;
p.marginright=5;
p.fontsize=10;
p.de.margin=4;
p(1).margintop=10;
p(2).margintop=6;

% p.fontsize=12;
p.select('all');
p(1,1).select();
ylabel({'Normalized cumulative','energy flux','[%]'});
set(gca,'XTick',10.^[3,4,5,5.477,6],'XLim',10.^[3 6],...
    'XTickLabel','',...
    'YLim',[0.1 110],'YTick',[ 1 10 20 30 40 50 60 70 80 90 100],...
    'YScale','linear','XScale','log');
title('PFISR Measurement');
grid on;
p(2,1).select();
ylabel({'Normalized cumulative','energy flux','[%]'});
set(gca,'XTick',10.^[3,4,5,5.477,6],'XLim',10.^[3 6],...
    'XTickLabel',{'1','10','100','300','1000'},...
    'YLim',[0.1 110],'YTick',[ 1 10 20 30 40 50 60 70 80 90 100],...
    'YScale','linear','XScale','log');
title('THEMIS Measurement');
grid on;
xlabel('Energy [keV]');
grid on;
% Plotting 1-D Cumulative Energy Flux

p(1,1).select();
plot(dataAvg.energyBin,...
    100*beforeOnset.dataAvg.cumuEnergyFlux./afterOnset.dataAvg.cumuEnergyFlux(end),...
     'k','LineWidth',1);
hold on;
medianValue=0.5*max(100*beforeOnset.dataAvg.cumuEnergyFlux./afterOnset.dataAvg.cumuEnergyFlux(end));
plot(10.^[3 6],[medianValue medianValue],'--k');
hold on;
plot(dataAvg.energyBin,...
    100*afterOnset.dataAvg.cumuEnergyFlux./afterOnset.dataAvg.cumuEnergyFlux(end),...
     'k','LineWidth',2);
 medianValue=0.5*max(100*afterOnset.dataAvg.cumuEnergyFlux./afterOnset.dataAvg.cumuEnergyFlux(end));
plot(10.^[3 6],[medianValue medianValue],'--k','LineWidth',2);
hold on;

legend('Before onset','Before onset median','After onset','After onset median','Location','southeast');
grid on;

p(2,1).select();

plot(thm.energyBin,...
    100*beforeOnset.thm.cumuEnergyFlux./afterOnset.thm.cumuEnergyFlux(end),...
     'r','LineWidth',1);
hold on;
medianValue=0.5*max(100*beforeOnset.thm.cumuEnergyFlux./afterOnset.thm.cumuEnergyFlux(end));
plot(10.^[3 6],[medianValue medianValue],'--r');
hold on;
plot(thm.energyBin,...
    100*afterOnset.thm.cumuEnergyFlux./afterOnset.thm.cumuEnergyFlux(end),...
     'r','LineWidth',2);
  medianValue=0.5*max(100*afterOnset.thm.cumuEnergyFlux./afterOnset.thm.cumuEnergyFlux(end));
plot(10.^[3 6],[medianValue medianValue],'--r','LineWidth',2);
hold on;
legend('Before onset','Before onset median','After onset','After onset median','Location','southeast');

grid on;

title('THEMIS-D Measurement');

%%
export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/cumlativeEnergyFlux.pdf' -pdf 
save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/cumlativeEnergyFlux.svg');
export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/cumlativeEnergyFlux.png' -r600 -png -nocrop
