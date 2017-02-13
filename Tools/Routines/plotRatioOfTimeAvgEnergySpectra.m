%%  Routine to Plot Cumulative Plots Energy Flux averaged across time
% We use SIC model, and constant alpha model to calculate q
% And the maximum entrop method to invert the spectra

clear all;
initialize_geodata;

%Global Variables
computer=getenv('computername');
if computer=='NITHIN-SURFACE'
    pfisrDataFileNameStr =...
    'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\PFISR_Energy_Spectra\Data\DataFile_2008_1.h5';
    thmDataFileNameStr =...
    'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Data_Mar_08_Event\space\Espectra_thm.mat';
    sicDataFileNameStr='C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Data\Ground\summaryresQinversion.mat';
else
    pfisrDataFileNameStr =...
    '/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/DataFile_2008_1.h5';
    thmDataFileNameStr =...
    '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Espectra_thm.mat';
    sicDataFileNameStr='/home/nithin/Documents/git-repos/ISR_School_2016/Antti  SIC Model/summaryresQinversion.mat';
end

timeMin=datenum('26 Mar 2008 08:00');
timeMax=datenum('26 Mar 2008 13:00');
altMin = 60;
altMax = 150;


load (thmDataFileNameStr);
% data_thm -> E, ebin, time
thm.energyFlux = data_thm.E(:,19:1:end)'*10^4; %[eV m-2 sr-1 s-1 eV-1]
thm.energyBin  = data_thm.ebin(19:1:end)';
thm.time       = data_thm.time;
clear data_thm;
energyBin = thm.energyBin;

%%  Mode 1: From measurements averaged across the PFISR beams
[ePopDenAvg, altPop, timePop] = get_pfisr_variable(pfisrDataFileNameStr, 'popl', 0);

% Tidying up Electron Density
[NeAvg, time] = crop_time(ePopDenAvg, timePop, timeMin, timeMax);
[NeAvg, alt] = crop_altitude(NeAvg, altPop,altMin, altMax);
NeAvg = interp_nans(NeAvg);
NeAvg(NeAvg<0)=10^9;

% Generate Production Rate
[qAvg,qAvgTime, alpha] = get_production_rate(NeAvg, alt, time, 2);

% Invert production rates
[dataAvg] = get_inverted_flux(qAvg, qAvgTime, alt, energyBin);

%%  Mode 1.1: Error in PFISR Measurements
[dePopDenAvg, altPop, timePop] = get_pfisr_variable(pfisrDataFileNameStr, 'dpopl', 0);

% Tidying up Electron Density
[dNeAvg, time] = crop_time(dePopDenAvg, timePop, timeMin, timeMax);
[dNeAvg, alt] = crop_altitude(dNeAvg, altPop,altMin, altMax);
dNeAvg = interp_nans(dNeAvg);
dNeAvg(dNeAvg<0)=10^9;

% Generate Production Rate
[dqAvg,dqAvgTime, alpha] = get_production_rate(dNeAvg, alt, time, 2);

% Invert production rates
[dataAvgError] = get_inverted_flux(dqAvg, dqAvgTime, alt, energyBin);
%% Mode 2: From SIC Measurements

load (sicDataFileNameStr);
[ionrate_raw, sic.alt] = crop_altitude(ionrate_raw, h,altMin, altMax);
[ionrate_smoothed, sic.alt] = crop_altitude(ionrate_smoothed, h,altMin, altMax);
sic.qRaw = ionrate_raw*10^6;
sic.qSmooth = ionrate_smoothed*10^6;
sic.Ne = Nedata*10^6;
sic.NeSIC = NeSIC*10^6;
sic.time = t;
clear h; clear ionrate_raw; clear ionrate_smoothed; clear Nedata; clear NeSIC; clear t; clear UT;

% Invert production rates
[dataAvgSIC] = get_inverted_flux(sic.qSmooth, sic.time, sic.alt, energyBin);

%% Calculating cumulative energy flux

beforeOnset.timeMin = datenum('26-Mar-2008 11:00');
beforeOnset.timeMax = datenum('26-Mar-2008 11:44');

beforeOnset.dataAvg.energyFlux = ...
    get_time_avg_time_series_data(dataAvg.energyFlux, dataAvg.time,...
    beforeOnset.timeMin, beforeOnset.timeMax);
beforeOnset.dataAvgError.energyFlux = ...
    get_time_avg_time_series_data(dataAvgError.energyFlux,...
    dataAvgError.time,beforeOnset.timeMin, beforeOnset.timeMax);
beforeOnset.dataAvgSIC.energyFlux = ...
    get_time_avg_time_series_data(dataAvgSIC.energyFlux,...
    dataAvgSIC.time,beforeOnset.timeMin, beforeOnset.timeMax);
beforeOnset.thm.energyFlux = ...
    get_time_avg_time_series_data(thm.energyFlux,...
    thm.time,beforeOnset.timeMin, beforeOnset.timeMax);

afterOnset.timeMin = datenum('26-Mar-2008 11:44');
afterOnset.timeMax = datenum('26-Mar-2008 12:00');

afterOnset.dataAvg.energyFlux = ...
    get_time_avg_time_series_data(dataAvg.energyFlux, dataAvg.time,...
    afterOnset.timeMin, afterOnset.timeMax);
afterOnset.dataAvgError.energyFlux = ...
    get_time_avg_time_series_data(dataAvgError.energyFlux,...
    dataAvgError.time,afterOnset.timeMin, afterOnset.timeMax);
afterOnset.dataAvgSIC.energyFlux = ...
    get_time_avg_time_series_data(dataAvgSIC.energyFlux,...
    dataAvgSIC.time,afterOnset.timeMin, afterOnset.timeMax);
afterOnset.thm.energyFlux = ...
    get_time_avg_time_series_data(thm.energyFlux,...
    thm.time,afterOnset.timeMin, afterOnset.timeMax);

%% Creating Panels
hFig=figure(3);
resize_figure(hFig);

clf
p=panel(hFig);
% p.pack(1,2);
panelSize = 60; %in mm
p.pack({{panelSize} {panelSize} {panelSize}},{{80}});
p.marginleft=25;
p.marginright=5;
p.fontsize=10;
p.de.margin=4;
p(1).margintop=10;
p(2).margintop=6;
p(3).margintop=6;
% p.fontsize=12;
p.select('all');
p(1,1).select();
ylabel({'\hspace{4pt} Diff. Energy Flux','[$eV m^{-2} sr^{-1} s^{-1} eV^{-1}$]'},'Interpreter','latex');
set(gca,'XTick',10.^[3,4,5,5.477,6],'XLim',10.^[3 6],...
    'XTickLabel','',...
    'YTick',10.^[6 7 8 9 10 11 12],'YLim',10.^[6 12],...
    'XTick',10.^[2,3,4,5,5.477,6],'YScale','log','XScale','log');
title('PFISR Measurement');
grid on;
p(2,1).select();
ylabel({'\hspace{4pt} Diff. Energy Flux','[$eV m^{-2} sr^{-1} s^{-1} eV^{-1}$]'},'Interpreter','latex');
set(gca,'XTick',10.^[3,4,5,5.477,6],'XLim',10.^[3 6],...
    'XTickLabel','',...
    'YTick',10.^[6 7 8 9 10 11 12],'YLim',10.^[6 12],...
    'XTick',10.^[2,3,4,5,5.477,6],'YScale','log','XScale','log');
title('THEMIS Measurement');
grid on;
p(3,1).select();
ylabel({'$$\int_{11:44}^{12:00} j(E) dt$$ / $$\int_{11:00}^{11:44} j(E) dt$$'},'Interpreter','latex');
set(gca, 'XTick',10.^[3,4,5,5.477,6],...
    'XTickLabel',{'1','10','100','300','1000'},...
    'XLim',10.^[3 6],...
    'YTick',10.^[-1 0 1 2 3],'YLim',10.^[-1 3],...
    'YTickLabel',{'0.1','1','10','100','1000'},...
    'YScale','log','XScale','log');
xlabel('Energy [keV]');
title('Ratio of spectra before and after onset');
grid on;

% Plotting Ratio of  Time Avg  Energy  Spectra

ratioPfisrEnergyFlux = afterOnset.dataAvg.energyFlux./beforeOnset.dataAvg.energyFlux;
ratioThmEnergyFlux = afterOnset.thm.energyFlux./beforeOnset.thm.energyFlux;

p(3,1).select();
plot(dataAvg.energyBin, ratioPfisrEnergyFlux,'--k');
hold on;
plot(thm.energyBin, ratioThmEnergyFlux,'-.r');
hold on;
text(-0.175,0.975,[char(96+3),')'],'Units','normalized','FontWeight','bold');

legend('Ionosphere (PFISR)', 'Plasmasheet (THEMIS)','Location','southeast');

% Plotting Time Avg  Energy  Spectra


p(1,1).select();
plot(dataAvg.energyBin, beforeOnset.dataAvg.energyFlux,'k','LineWidth',1);
hold on;
plot(dataAvg.energyBin, afterOnset.dataAvg.energyFlux,'k','LineWidth',2);
hold on;
text(-0.175,0.975,[char(96+1),')'],'Units','normalized','FontWeight','bold');

legend(char({'Before Onset'}),...
    char({'After Onset'}),'Location', 'southwest');

p(2,1).select();
plot(thm.energyBin, beforeOnset.thm.energyFlux,'r','LineWidth',1);
hold on;
plot(thm.energyBin, afterOnset.thm.energyFlux,'r','LineWidth',2);
hold on;
text(-0.175,0.975,[char(96+2),')'],'Units','normalized','FontWeight','bold');

legend(char({'Before Onset'}),...
    char({'After Onset'}),'Location', 'southwest');

%% Plot
if computer=='NITHIN-SURFACE'
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\RatioBeforeAfterOnset.pdf' -pdf -nocrop
    save('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\RatioBeforeAfterOnset.svg');
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\RatioBeforeAfterOnset.png' -r600 -png -nocrop
else
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/RatioBeforeAfterOnset.pdf' -pdf 
    save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/RatioBeforeAfterOnset.svg');
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/RatioBeforeAfterOnset.png' -r600 -png -nocrop
end
