%% Inputs
tic
% PFISR

% pfisrExpFileName = '20080215.005_bc_2min-Ne-cal.h5';
pfisrExpFileName = '20080326.001_bc_2min-fitcal.h5';

if isunix
    baseDir = '/media/nithin/PFISR_002_006/';
    pfisrRootPath = [baseDir,'PFISR Processed/Sporadic04/'];
    dascRootPath = [baseDir,'DASC/'];
    outputH5FileStr = [baseDir,'PFISR Processed/Sporadic04/',...
    pfisrExpFileName(1:20),'-energyFlux.h5'];
    outputFiguresFolder = ['/media/nithin/PFISR_002_006/PFISR Processed/Sporadic04/Figures_',pfisrExpFileName(1:8),'/'];
else
    baseDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Event 1\';
    pfisrRootPath = baseDir;
    dascRootPath = 'C:\Users\nithin\Documents\GitHub\LargeFiles\DASC\';
    outputH5FileStr = [baseDir,...
        pfisrExpFileName(1:20),'-energyFlux_v2.h5'];
    outputFiguresFolder = [baseDir,'Figures_v2_',pfisrExpFileName(1:8),'\'];
end

% minTimeStr = [];
% maxTimeStr = [];

minTimeStr = '26-Mar-2008 11:00';
maxTimeStr = '26-Mar-2008 12:00';

pfisrFileNameStr = [pfisrRootPath,pfisrExpFileName];
minAlt = 60;
maxAlt = 200;
projectionAlt = 60; % km
nEnergyBins = 30;
minE = 10^3;
maxE = 10^6;

% DASC
dascMinElevation = 0;
dascCalFileAz = [];
dascCalFileEl = [];
dascSetDownloadFlag = true;

%Images
energySlice = 3; %keV
energyFluxLim = [10 13];
timeMinStr = [];
timeMaxStr = [];
latLim = [];
lonLim  = [];
setStoreImage = true;

% Files
% Plotting/Videos
outputVideoStr = [pfisrExpFileName(1:8),'_energyFlux_',num2str(energySlice),'keV_2min.avi'];


%% Creating H5 File
multiWaitbar('Close All');

if ~isfile(outputH5FileStr)
    [status]=generate_energy_spectra_data_product(pfisrExpFileName,pfisrRootPath,...
        dascRootPath,projectionAlt,logspace(log10(minE),log10(maxE),nEnergyBins),...
        [minAlt maxAlt], outputH5FileStr, minTimeStr, maxTimeStr,...
        dascMinElevation,dascCalFileAz,dascCalFileEl, dascSetDownloadFlag);
end
set(0, 'DefaultFigureVisible', 'off');
create_energy_spectra_images(outputH5FileStr,...
    outputFiguresFolder,outputVideoStr,energySlice,energyFluxLim,...
    minTimeStr,maxTimeStr,latLim,lonLim,setStoreImage);
set(0, 'DefaultFigureVisible', 'on');
toc