%% For Paper 2, Figure 2
% Time variation of EEA, through FYKN, pokerFlat, and gako
clear all;

%% Load THEMIS foot prints
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');

% Load omni  & magnetic field model
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 09:00','26 Mar 2008 12:00');
magFieldNo = find_irbem_magFieldModelNo('TS96');
maginput = filter_irbem_maginput(magFieldNo,maginput);

% Load PFISR data
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';

% Load POES/NOAA-17
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);

%% Plot figure
latLim = [61,68];
lonLim = [-153,-141];
deltaLat = 2;
deltaLon = 5;
storeImageDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Figures\Draft\Figure2\';

% Plotting
hFig1=figure('visible','on');
timeStr = '26 Mar 2008 10:21:50';
time = datenum(timeStr);
[ax1]=combine_2D_plots_v2(fileStr,hFig1,...
      'maps',{'OpticalImage','OpticalImage','OpticalImage','EnergyFluxMap'},...
      'sites',{'pokerFlat','gako','fykn','pokerFlat'},...
      'thisTime',time,...
      'latLim',latLim,...
      'lonLim',lonLim,...
      'elCutOff',30,...
      'deltaLat',deltaLat,...
      'deltaLon',deltaLon,...
      'opticalLim',[250 400],...
      'peakIonizationAltitude',85,...
      'setStoreImage',false);
% resize_figure(hFig1,279.4/3,215.9/3);

thdTimeIndx = find_time(padData.thd.time,timeStr);
theTimeIndx = find_time(padData.the.time,timeStr);
hold on;
scatterm(padData.thd.latFoot(thdTimeIndx), padData.thd.lonFoot(thdTimeIndx),...
    '^b','filled');
scatterm(padData.the.latFoot(theTimeIndx), padData.the.lonFoot(theTimeIndx),...
    'sc','filled');
%%
hFig2=figure('visible','on');
% resize_figure(hFig2,279.4/3,215.9/3);
timeStr = '26 Mar 2008 11:18:50';
time = datenum(timeStr);
[ax2]=combine_2D_plots_v2(fileStr,hFig2,...
      'maps',{'OpticalImage','OpticalImage','OpticalImage','EnergyFluxMap'},...
      'sites',{'gako','fykn','pokerFlat','pokerFlat'},...
      'thisTime',time,...
      'latLim',latLim,...
      'lonLim',lonLim,...
      'elCutOff',30,...
      'deltaLat',deltaLat,...
      'deltaLon',deltaLon,...
      'opticalLim',[250 400],...
      'peakIonizationAltitude',85,...
      'setStoreImage',false);

thdTimeIndx = find_time(padData.thd.time,timeStr);
theTimeIndx = find_time(padData.the.time,timeStr);


hFig3=figure('visible','on');
% resize_figure(hFig3,279.4/3,215.9/3);
timeStr = '26 Mar 2008 11:35:50';
time = datenum(timeStr);
[ax3]=combine_2D_plots_v2(fileStr,hFig3,...
      'maps',{'OpticalImage','OpticalImage','OpticalImage','EnergyFluxMap'},...
      'sites',{'fykn','pokerFlat','gako','pokerFlat'},...
      'thisTime',time,...
      'latLim',latLim,...
      'lonLim',lonLim,...
      'elCutOff',30,...
      'deltaLat',deltaLat,...
      'deltaLon',deltaLon,...
      'opticalLim',[250 400],...
      'peakIonizationAltitude',85,...
      'setStoreImage',false);
