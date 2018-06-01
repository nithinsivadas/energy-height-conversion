%% Inputs
tic
% PFISR

% pfisrExpFileName = '20080215.005_bc_2min-Ne-cal.h5';
pfisrExpFileName = '20080326.001_bc_2min-fitcal.h5';
pfisrdTime = pfisrExpFileName(regexp(pfisrExpFileName,'bc_')+3:regexp(pfisrExpFileName,'-fit')-1);
if isunix
    baseDir = '/media/nithin/PFISR_002_006/';
    pfisrRootPath = [baseDir,'PFISR Processed/Sporadic04/'];
    dascRootPath = [baseDir,'DASC/'];
    outputH5FileStr = [baseDir,'PFISR Processed/Sporadic04/',...
    pfisrExpFileName(1:regexp(pfisrExpFileName,'-fit')),'energyFlux_v7_1.h5'];
    outputFigureBaseDir = '/media/nithin/PFISR_002_006/PFISR Processed/Sporadic04/';
    outputFigureFolderStr = ['Figures_v7_',pfisrExpFileName(1:8)];
else
    baseDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Event 1\';
    pfisrRootPath = baseDir;
    dascRootPath = 'C:\Users\nithin\Documents\GitHub\LargeFiles\DASC\';
    outputH5FileStr = [baseDir,...
        pfisrExpFileName(1:20),'-energyFlux_v1.h5'];
    outputFiguresFolder = [baseDir,'Figures_v1_',pfisrExpFileName(1:8),'\'];
end

% minTimeStr = [];
% maxTimeStr = [];

minTimeStr = '26 Mar 2008 11:20';
maxTimeStr = '26 Mar 2008 11:30';

pfisrFileNameStr = [pfisrRootPath,pfisrExpFileName];
minAlt = 60;
maxAlt = 200;
projectionAltPFISR = 60; % km, altitude of origin of magnetic field aligned lines
nEnergyBins = 30;
minE = 10^3;
maxE = 10^6;

% DASC
dascMinElevation = 0;
dascCalFileAz = [];
dascCalFileEl = [];
dascSetDownloadFlag = true;
projectionAltDASC = 110; %110 km

%Images
% energySlice = [3, 5, 10, 30, 100]; %keV
% energyFluxLim = [10 13;9 11; 9 11; 8 10; 8 10];
energySlice = 3;
energyFluxLim = [10 13];
opticaLim = [300 600];
timeMinStr = [];
timeMaxStr = [];
latLim = [];
lonLim  = [];
setStoreImage = true;




%% Creating H5 File
multiWaitbar('Close All');

if ~isfile(outputH5FileStr)
    [status]=generate_energy_spectra_data_product(pfisrExpFileName,pfisrRootPath,...
        dascRootPath,projectionAltPFISR,projectionAltDASC,logspace(log10(minE),log10(maxE),nEnergyBins),...
        [minAlt maxAlt], outputH5FileStr, minTimeStr, maxTimeStr,...
        dascMinElevation,dascCalFileAz,dascCalFileEl, dascSetDownloadFlag);
end
set(0, 'DefaultFigureVisible', 'off');
multiWaitbar('Energy slice progress',0);
nEnergy = length(energySlice);
for iEnergy = 1:1:nEnergy 
% Plotting/Videos
    outputVideoStr = [pfisrExpFileName(1:8),'_energyFlux_',num2str(energySlice(iEnergy)),'keV_',pfisrdTime,'.avi'];
    outputFiguresFolder = [outputFigureBaseDir,outputFigureFolderStr,'_',num2str(energySlice(iEnergy)),'keV_',pfisrdTime,'/'];
    create_energy_spectra_images(outputH5FileStr,...
        outputFiguresFolder,outputVideoStr,energySlice(iEnergy),energyFluxLim(iEnergy,:),opticalLim,...
        minTimeStr,maxTimeStr,latLim,lonLim,setStoreImage);
    multiWaitbar('Energy slice progress','Increment',1./nEnergy);
end

set(0, 'DefaultFigureVisible', 'on');
toc