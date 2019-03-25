%% Create HDF5 file containing
% PFISR
% DASC, THEMIS ASIs
% THEMIS S/C & NOAA Lat, Lon Values
% Conductivity

%% Initialization
% pfisrExpFileName = '20080326.001_bc_15sec-fitcal.h5';
pfisrExpFileName = '20101018.001_bc_15sec-fitcal.h5';
pfisrdTime = pfisrExpFileName(regexp(pfisrExpFileName,'bc_')+3:regexp(pfisrExpFileName,'-fitcal')-1);

baseDir = '/media/nithin/PFISR_002_006/';
base2Dir = 'PFISR Processed/Event_List/';

pfisrRootPath = [baseDir,base2Dir];
dascRootPath = [baseDir,'DASC/'];
outputH5FileStr = [baseDir,base2Dir,...
pfisrExpFileName(1:regexp(pfisrExpFileName,'-fitcal')),'full_v1.h5'];
pfisrFileNameStr = [pfisrRootPath,pfisrExpFileName];

% Time
minTimeStr = '18 Oct 2010 02:00';
maxTimeStr = '18 Oct 2010 18:00';
dateFormat = 'dd mmm yyyy HH:SS';
% Altitude
minAlt = 60;
maxAlt = 200;
projectionAltPFISR = 60; % km, altitude of origin of magnetic field aligned lines

% Energy
nEnergyBins = 30;
minE = 10^3;
maxE = 10^6;

% DASC
dascMinElevation = 0;
dascCalFileAz = [];
dascCalFileEl = [];
dascSetDownloadFlag = true;
projectionAltDASC = 110; %110 km

% NOAA Mat File
% noaaMatFileStr = '/media/nithin/PFISR_002_006/Nithin/NOAA17/n1720080326.mat';
omniH5FileStr = '/home/nithin/Documents/git-repos/LargeFiles/omni/omni.h5';

%% Write PFISR Energy Spectra, and DASC data
[status]=generate_energy_spectra_data_product(pfisrExpFileName,pfisrRootPath,...
        dascRootPath,projectionAltPFISR,projectionAltDASC,logspace(log10(minE),log10(maxE),nEnergyBins),...
        [minAlt maxAlt], outputH5FileStr, minTimeStr, maxTimeStr,...
        dascMinElevation,dascCalFileAz,dascCalFileEl, dascSetDownloadFlag);
fprintf('PFISR Energy & DASC Images: Done \n');

%% Write THM ASI
% sites={'gako','fykn','inuv','whit','mcgr','kian'};
sites = {'mcgr','kian'};
multiWaitbar('Processing different cameras',0);
id = 1./length(sites);
for i=1:1:length(sites)
multiWaitbar('Processing different cameras','Increment',id);
siteName = sites{i};
disp(['Processing ',upper(siteName)]);
status = create_thg_hdf5(siteName,outputH5FileStr);
end
fprintf('THM ASIs: Done \n');

%% [Write a better version, that has a better chunkSize and dataSize definition]
% Write Magnetic Field Model 
% magFieldModelStr = {'TS89','TS96','TS01'};
% id = 1./length(magFieldModelStr);
% for i = 1:1:length(magFieldModelStr)
%     magFieldModelNo = find_irbem_magFieldModelNo(magFieldModelStr{i});
%     [status] = create_magneticFieldProjection_hdf5(magFieldModelNo,...
%         outputH5FileStr,omniH5FileStr,'options',[0,0,0,0,0],'pixels',32);
% end

% fprintf('Tsyganenko Models: Done \n');
%% Write thm spacecraft data to hdf5
sc = {'tha','thb','thc','thd','the'};
id = 1./length(sc);
for i = 1:1:length(sc)
    add_thm_hdf5(sc{i},outputH5FileStr,omniH5FileStr,...
        'minTimeStr',minTimeStr,...
        'maxTimeStr',maxTimeStr,...
        'dateFormat',dateFormat,...
        'magneticFieldModel','TS89');
end

fprintf('THEMIS State Data: Done \n');
%% Write noaa17
% add_temp_noaa17_hdf5(noaaMatFileStr, outputH5FileStr, omniH5FileStr,...
%         'minTimeStr',minTimeStr,...
%         'maxTimeStr',maxTimeStr,...
%         'dateFormat',dateFormat,...
%         'magneticFieldModel','TS89'); 
%     
% fprintf('NOAA17 State Data: Done \n');
%% Write conductivity
add_conductivity_hdf5(outputH5FileStr,outputH5FileStr,true);

fprintf('Conductivity Estimates: Done \n');