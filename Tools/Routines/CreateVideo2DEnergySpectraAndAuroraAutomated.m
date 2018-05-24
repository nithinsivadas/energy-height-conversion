%% Inputs
tic
% PFISR

pfisrExpFileName = '20080215.005_bc_2min-Ne-cal.h5';
% pfisrExpFileName = '20080326.001_bc_15sec-fitcal.h5';

if isunix
    baseDir = '/media/nithin/PFISR_002_006/';
    pfisrRootPath = [baseDir,'PFISR Processed/Sporadic03/'];
    dascRootPath = [baseDir,'DASC'];
    outputH5FileStr = [baseDir,'PFISR Processed/Sporadic03/',...
    pfisrExpFileName(1:20),'-energyFlux.h5'];
    outputFiguresFolder = ['/media/nithin/PFISR_002_006/PFISR Processed/Sporadic03/Figures_',pfisrExpFileName(1:8),'/'];
else
    baseDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Event 1\';
    pfisrRootPath = baseDir;
    dascRootPath = 'C:\Users\nithin\Documents\GitHub\LargeFiles\DASC';
    outputH5FileStr = [baseDir,...
        pfisrExpFileName(1:20),'-energyFlux.h5'];
    outputFiguresFolder = [baseDir,'Figures_',pfisrExpFileName(1:8),'\'];
end

minTimeStr = [];
maxTimeStr = [];
pfisrFileNameStr = [pfisrRootPath,pfisrExpFileName];
minAlt = 60;
maxAlt = 200;
projectionAlt = 70; % km
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
energyFluxLim = [9 12];
timeMinStr = [];
timeMaxStr = [];
latLim = [];
lonLim  = [];
setStoreImage = true;

% Files
% Plotting/Videos
outputVideoStr = [pfisrExpFileName(1:8),'_energyFlux_',num2str(energySlice),'keV.avi'];


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