%%  Routine to Plot PFISR based energy spectra in varous modes

clear all;
initialize_geodata;

%Global Variables
computer=getenv('computername');
if computer=='NITHIN-SURFACE'
   pfisrDataFileNameStr =...
    'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\PFISR_Energy_Spectra\Data\DataFile_2008_1.h5';
   thmDataFileNameStr =...
    'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Data_Mar_08_Event\space\Espectra_thm.mat';
else
   pfisrDataFileNameStr =...
    '/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/DataFile_2008_1.h5';
   thmDataFileNameStr =...
    '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Espectra_thm.mat';
end

timeMin=datenum('26 Mar 2008 08:00');
timeMax=datenum('26 Mar 2008 13:00');
altMin = 60;
altMax = 180;


load (thmDataFileNameStr);
% data_thm -> E, ebin, time
thm.energyFlux = data_thm.E(:,19:1:end)'*10^4; %[eV m-2 sr-1 s-1 eV-1]
thm.energyBin  = data_thm.ebin(19:1:end)';
thm.time       = data_thm.time;
clear data_thm;
energyBin = thm.energyBin;

%% Mode 1: Inversion from measurements along the magnetic field line
[ePopDenMag, altPop, timePop] = get_pfisr_variable(pfisrDataFileNameStr, 'popl', -1);

% Tidying up Electron Density
[NeMag, time] = crop_time(ePopDenMag, timePop, timeMin, timeMax);
[NeMag, alt] = crop_altitude(NeMag, altPop,altMin, altMax);
NeMag = interp_nans(NeMag);
NeMag(NeMag<0)=10^9;

% Generate Production Rate
[qMag,qMagTime, alpha] = get_production_rate(NeMag, alt, time, 2);

% Invert production rates
[dataMag] = get_inverted_flux(qMag, qMagTime, alt, energyBin);
dataMag.energyFlux(dataMag.energyFlux<0)=10^6;


%%  Mode 2: From measurements averaged across the PFISR beams
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


%% Plotting 2-D Electron Density neasurements
figure;

subplot (2,1,1)
plot_2D_time_series(time, alt, log10(NeMag), 0.5, 1);
title('Uncorrected Electron Density: Along B Field alinged Beam 1 log_1_0 [m^-^3]');
caxis([6, 12]);

subplot (2,1,2)
plot_2D_time_series(time, alt, log10(NeAvg), 0.5, 1);
title('Uncorrected Electron Density: Averaged across all beams log_1_0 [m^-^3]');
caxis([6, 12]);

%% Plotting 2-D Energy Spectra measurements
figure;
subplot (3,1,1)
plot_2D_time_series(dataMag.time, energyBin, log10(dataMag.energyFlux), 0.5, 3);
title('Differential Energy Flux: Along B Field alinged Beam 1 log_1_0 [eV m^-^2 sr^-^1 s^-^1 eV^-^1]');
caxis([6, 12]);

subplot (3,1,2)
plot_2D_time_series(dataAvg.time, energyBin, log10(dataAvg.energyFlux), 0.5, 3);
title('Differential Energy Flux: Averaged across all beams log_1_0 [eV m^-^2 sr^-^1 s^-^1 eV^-^1]');
caxis([6, 12]);

subplot (3,1,3)
plot_2D_time_series(thm.time, energyBin, log10(thm.energyFlux), 0.5, 3);
title('Differential Energy Flux: THEMIS D Spacecraft log_1_0 [eV m^-^2 sr^-^1 s^-^1 eV^-^1]');
caxis([6, 12]);

%% Plotting 1-D Validation Plots
thm.alt = dataAvg.alt;
thm.qForward = dataAvg.A*energy_to_num(thm.energyFlux, thm.time, thm.energyBin);
thm.NeForward = q_to_Ne(thm.qForward, thm.alt, thm.time);

%%
figure;
subplot(1,3,1)
p(1)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, NeAvg, '26-Mar-2008 09:47', 1);
p(1).Color = 'b';
hold on;
% p(2)=plot_1D_time_slice(dataMag.time, dataMag.alt, NeMag, '26-Mar-2008 09:47', 1);
% p(2).Color = 'r'; p(2).LineStyle = '--';p(2).LineWidth = 1;
% hold on;
p(3)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, q_to_Ne(dataAvg.qInverted, dataAvg.alt, dataAvg.time), '26-Mar-2008 09:47', 1);
p(3).Color = 'k'; p(3).LineStyle = '-.';p(3).LineWidth = 1;
hold on;
p(4)=plot_1D_time_slice(thm.time, thm.alt, thm.NeForward, '26-Mar-2008 09:47', 1);
p(4).Color = 'g'; p(4).LineStyle = ':';p(4).LineWidth = 2;
xlim([10^8 10^12]);
legend('N_e_ _P_F_I_S_R Measured','N_e_ _P_F_I_S_R Inverted', 'N_e_ _T_H_E_M_I_S');

subplot(1,3,2)
p(1)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, qAvg, '26-Mar-2008 09:47', 2);
p(1).Color = 'b';
hold on;
% p(2)=plot_1D_time_slice(dataMag.time, dataMag.alt, qMag, '26-Mar-2008 09:47', 2);
% p(2).Color = 'r'; p(2).LineStyle = '--';p(2).LineWidth = 1;
% hold on;
p(3)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, dataAvg.qInverted, '26-Mar-2008 09:47', 2);
p(3).Color = 'k'; p(3).LineStyle = '-.';p(3).LineWidth = 1;
p(4)=plot_1D_time_slice(thm.time, thm.alt, thm.qForward, '26-Mar-2008 09:47',2);
p(4).Color = 'g'; p(4).LineStyle = ':';p(4).LineWidth = 2;
xlim([10^5 10^11]);
legend('q_P_F_I_S_R Measured','q_P_F_I_S_R Inverted', 'q_T_H_E_M_I_S');
title(datestr(dataAvg.time(find_time(dataAvg.time, ('26-Mar-2008 09:47')))));


subplot(1,3,3)
p(1)=plot_1D_time_slice(dataAvg.time, dataAvg.energyBin, dataAvg.energyFlux, '26-Mar-2008 09:47', 3);
p(1).Color = 'b';
hold on;
% p(2)=plot_1D_time_slice(dataMag.time, dataMag.energyBin, dataMag.energyFlux, '26-Mar-2008 09:47', 3);
% p(2).Color = 'r'; p(2).LineStyle = '--';p(2).LineWidth = 1;
% hold on;
p(3)=plot_1D_time_slice(thm.time, thm.energyBin, thm.energyFlux, '26-Mar-2008 09:47', 3);
p(3).Color = 'g'; p(3).LineStyle = ':';p(3).LineWidth = 2;
% xlim([10^5 10^11]);
legend('j(E)_P_F_I_S_R ','j(E)_T_H_E_M_I_S');