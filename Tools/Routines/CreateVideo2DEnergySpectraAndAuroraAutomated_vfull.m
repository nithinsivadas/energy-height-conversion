%% Inputs

% PFISR
pfisrExpFileName = '20080215.005_bc_2min-Ne-cal.h5';
minTimeStr = [];
maxTimeStr = [];
pfisrRootPath = '/media/nithin/PFISR_002_006/PFISR Processed/Sporadic03/';
pfisrFileNameStr = [pfisrRootPath,pfisrExpFileName];
minAlt = 60;
maxAlt = 200;
projectionAlt = 70; % km
nEnergyBins = 30;
minE = 10^3;
maxE = 10^6;

% DASC
dascRootPath = '/media/nithin/PFISR_002_006/DASC';
dascMinElevation = 30;
dascCalFileAz = [];
dascCalFileEl = [];
dascSetDownloadFlag = false;


outputH5FileStr = ['/media/nithin/PFISR_002_006/PFISR Processed/Sporadic03/',...
    pfisrExpFileName(1:20),'-energyFlux_full.h5'];

%% Creating H5 File
[data,status]=generate_energy_spectra_data_product(pfisrExpFileName,pfisrRootPath,...
    dascRootPath,projectionAlt,logspace(log10(minE),log10(maxE),nEnergyBins),...
    [minAlt maxAlt], outputH5FileStr, minTimeStr, maxTimeStr,...
    dascMinElevation,dascCalFileAz,dascCalFileEl, dascSetDownloadFlag);
