%% Inputs
tic
% PFISR

% pfisrExpFileName = '20080215.005_bc_2min-Ne-cal.h5';
% pfisrExpFileName = '20080326.001_bc_2min-fitcal.h5';
% pfisrExpFileName = '20100528.001_bc_2min-Ne-cal.h5';
pfisrExpFileName = '20101018.001_bc_15sec-fitcal.h5';
% pfisrExpFileName = '20180911.001_bc_1min-fitcal.h5';
% pfisrExpFileName = '20080326.001_bc_15sec-fitcal.h5';
% pfisrdTime = pfisrExpFileName(regexp(pfisrExpFileName,'bc_')+3:regexp(pfisrExpFileName,'-Ne')-1);
pfisrdTime = pfisrExpFileName(regexp(pfisrExpFileName,'bc_')+3:regexp(pfisrExpFileName,'-fitcal')-1);

baseDir = '/media/nithin/PFISR_002_006/';
base2Dir = 'PFISR Processed/Event_List/';
% base2Dir = 'PFISR Processed/ThemisD1.v01/20180911.001/';
pfisrRootPath = [baseDir,base2Dir];
dascRootPath = [baseDir,'DASC/'];
outputH5FileStr = [baseDir,base2Dir,...
pfisrExpFileName(1:regexp(pfisrExpFileName,'-fitcal')),'full_v110.h5'];
% pfisrExpFileName(1:regexp(pfisrExpFileName,'-Ne')),'energyFlux.h5'];

outputFigureBaseDir = '/media/nithin/PFISR_002_006/PFISR Processed/Event_List/';
outputFigureFolderStr = ['Figures_',pfisrExpFileName(1:8)];


% minTimeStr = [];
% maxTimeStr = [];

minTimeStr = '18 Oct 2010 07:30';
maxTimeStr = '18 Oct 2010 10:30';

pfisrFileNameStr = [pfisrRootPath,pfisrExpFileName];
minAlt = 60;
maxAlt = 200;
% maxAlt = 120;
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
energySlice = 100;
energyFluxLim = [8 10];
opticalLim = [300 600];
timeMinStr = [];
timeMaxStr = [];
latLim = [];
lonLim  = [];
setStoreImage = true;




%% Creating H5 File
multiWaitbar('Close All');

% if ~isfile(outputH5FileStr)
    [status]=generate_energy_spectra_data_product(pfisrExpFileName,pfisrRootPath,...
        dascRootPath,projectionAltPFISR,projectionAltDASC,logspace(log10(minE),log10(maxE),nEnergyBins),...
        [minAlt maxAlt], outputH5FileStr, minTimeStr, maxTimeStr,...
        dascMinElevation,dascCalFileAz,dascCalFileEl, dascSetDownloadFlag);
% end

% sites = {'gako'};
% sites={'gako','fykn','inuv','whit','mcgr','kian'};
% multiWaitbar('Processing different cameras',0);
% for i=1:1:length(sites)
% multiWaitbar('Processing different cameras','Increment',1./length(sites));
% % outputFileStr ='/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-energyFlux_v85.h5';
% siteName = sites{i};
% disp(['Processing ',upper(siteName)]);
% status = create_thg_hdf5(siteName,outputFileStr);
% end
% multiWaitbar('Close All');

% set(0, 'DefaultFigureVisible', 'off');
% multiWaitbar('Energy slice progress',0);
% nEnergy = length(energySlice);
% for iEnergy = 1:1:nEnergy 
% % Plotting/Videos
%     outputVideoStr = [pfisrExpFileName(1:8),'_energyFlux_',num2str(energySlice(iEnergy)),'keV_',pfisrdTime,'.avi'];
%     outputFiguresFolder = [outputFigureBaseDir,outputFigureFolderStr,'_',num2str(energySlice(iEnergy)),'keV_',pfisrdTime,'/'];
%     create_energy_spectra_images(outputH5FileStr,...
%         outputFiguresFolder,outputVideoStr,energySlice(iEnergy),energyFluxLim(iEnergy,:),opticalLim,...
%         minTimeStr,maxTimeStr,latLim,lonLim,setStoreImage);
%     multiWaitbar('Energy slice progress','Increment',1./nEnergy);
% end
% 
% set(0, 'DefaultFigureVisible', 'on');
toc