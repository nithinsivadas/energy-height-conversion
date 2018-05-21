%% Inputs
tic
% PFISR
pfisrExpFileName = '20080215.005_bc_2min-Ne-cal.h5';
minTimeStr = '15 Feb 2008 15:40';
maxTimeStr = '16 Feb 2008 04:50';
% pfisrRootPath = '/media/nithin/PFISR_002_006/PFISR Processed/Sporadic03/';
pfisrRootPath = 'G:\My Drive\Research\Projects\Paper 2\Data\Madrigal Events\';
pfisrFileNameStr = [pfisrRootPath,pfisrExpFileName];
minAlt = 60;
maxAlt = 200;
projectionAlt = 70; % km
nEnergyBins = 30;
minE = 10^3;
maxE = 10^6;

% DASC
% dascRootPath = '/media/nithin/PFISR_002_006/DASC';
dascRootPath = 'C:\Users\nithin\Documents\GitHub\LargeFiles\DASC';
dascMinElevation = 0;
dascCalFileAz = [];
dascCalFileEl = [];
dascSetDownloadFlag = false;

% Files
outputH5FileStr = ['G:\My Drive\Research\Projects\Paper 2\Data\Madrigal Events\',...
    pfisrExpFileName(1:20),'-energyFlux_sample_3.h5'];
% outputH5FileStr = ['/media/nithin/PFISR_002_006/PFISR Processed/Sporadic03/',...
%     pfisrExpFileName(1:20),'-energyFlux_sample.h5'];
% Plotting/Videos
outputFiguresFolder = 'G:\My Drive\Research\Projects\Paper 2\Data\Madrigal Events\Figures\';
outputVideoStr = '20080216_energyFlux_3keV.avi';
energySlice = 3; %keV
energyFluxLim = [9 12];
timeMinStr = [];
timeMaxStr = [];
latLim = [];
lonLim  = [];
setStoreImage = true;

%% Creating H5 File
if ~isfile(outputH5FileStr)
    [data,status]=generate_energy_spectra_data_product(pfisrExpFileName,pfisrRootPath,...
        dascRootPath,projectionAlt,logspace(log10(minE),log10(maxE),nEnergyBins),...
        [minAlt maxAlt], outputH5FileStr, minTimeStr, maxTimeStr,...
        dascMinElevation,dascCalFileAz,dascCalFileEl, dascSetDownloadFlag);
end
% 
% [time,timePrimary] = create_energy_spectra_images(outputH5FileStr,...
%     outputFiguresFolder,outputVideoStr,energySlice,energyFluxLim,...
%     timeMinStr,timeMaxStr,latLim,lonLim,setStoreImage);
toc