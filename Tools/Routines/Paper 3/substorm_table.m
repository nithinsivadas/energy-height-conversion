%% Create a substorm table
% This routine creates a table of superMag substorms, and connects its time
% of arrival with the time ranges where data is available from PFISR and/or
% DASC. 
% Input:
%       amisrWebDatabase.h5 (PFISR database)
%       substorms_superMag_20201130.txt
%       omni.h5
% Output: 
%       table_of_substorms_full.ma
%       T - Table with all substorms which have PFISR data
%       T1 - Table with closest substorms that have PFSIR and DASC data

clear all;

if strcmp(get_computer_name,'nithin-carbon')
    dataDir = 'G:\My Drive\Research\Projects\Data\';
    storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\';
elseif strcmp(get_computer_name,'aurora1-optiplex-780')
    dataDir = '/media/nithin/Elements/Nithin/Data/';
    storeDir = '/media/nithin/Elements/Nithin/Data/Paper_3/';
elseif strcmp(get_computer_name,'scc-lite')
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/scratch/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
else
%     error(['Not configured for this computer: ',get_computer_name]);
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/projectnb/semetergrp/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
end

outputAMISRFileStr = 'amisrWebDatabase.h5';
amisrDatabaseStr = [dataDir,outputAMISRFileStr];

superMagFileStr = [storeDir,'substorms_superMag_20201130.txt'];

dascFileStr = [storeDir,'dascDatabase.h5'];
omniFileStr = [dataDir,'omni.h5'];

timeMinStr = "01 Dec 2006";
timeMaxStr = "31 Dec 2019";

% Substorms at PFISR [IMPORTANT] : Range of closeness to PFISR
Dmlt = 2;
Dmlat = 50; %desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat

% Poker Flat geolocation
pkrGLAT = 65.126;
pkrGLON = -147.47;
pkrh0=0.693;

tic 
[T,T1, superMag] = substorm_create_table(dataDir,outputAMISRFileStr,amisrDatabaseStr,...
  superMagFileStr, dascFileStr, omniFileStr,...
  timeMinStr, timeMaxStr, Dmlt, Dmlat, pkrGLAT, pkrGLON, pkrh0);
toc 

disp(['Parameters of the run: Dmlt = ',num2str(Dmlt),' Dmlat = ',num2str(Dmlat)]);
disp(['The .mat file is stored in: ',storeDir,'table_of_substorms_full.mat']);

save([storeDir,'table_of_substorms_full.mat'],'T1','T','superMag','Dmlat','Dmlt','pkrGLAT','pkrGLON','dascFileStr','amisrDatabaseStr','omniFileStr','superMagFileStr');

function [T, T1, superMag] = substorm_create_table(dataDir,outputAMISRFileStr,amisrDatabaseStr,...
  superMagFileStr, dascFileStr, omniFileStr,...
  timeMinStr, timeMaxStr, Dmlt, Dmlat, pkrGLAT, pkrGLON, pkrh0)
% Loading database
% Create amisr database
if ~isfile(amisrDatabaseStr)
    status = create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,[dataDir,outputAMISRFileStr]);
end
%
% Load supermag database
superMag = extract_superMag_data(superMagFileStr);
superMag.stormID = (1:length(superMag.time))';

% Load amisr data
amisr = extract_amisr_data(amisrDatabaseStr);

% Load omni data
omni.AE = h5read(omniFileStr,'/Indices/AE');
omni.time = unixtime2matlab(h5read(omniFileStr,'/Time'));

% Calculations
% Estimating PFISR magnetic coordinates 
[superMag.pfisrMlat,superMag.pfisrMlon,superMag.pfisrMlt] = get_magnetic_coordinates([pkrh0,pkrGLAT,pkrGLON],superMag.time(:));

% Interpolating AE index of the substorms
superMag.AE = interp1(omni.time,omni.AE,superMag.time);

%%
% Selecting substorms closest to PFISR location
deltaMLT = absDiffMLT(superMag.pfisrMlt,superMag.mlt);
desiredMLTIndx = abs(deltaMLT)<Dmlt;
desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat;
closestSubstormIndx = desiredMLTIndx & desiredMLATIndx; 

%% Adding the PFISR experiments running during the substorm time

expBCArray = barker_coded_experiments();
% Barker Coded PFISR Experiment Filter
bcFilterIndx = zeros(1,length(amisr.expId));
for iexp = 1:1:length(expBCArray)
    bcFilterIndx = bcFilterIndx|strcmp(strtrim(expBCArray(iexp)),cellstr(deblank(amisr.expName)));
end

for iStorm = 1:1:length(superMag.stormID)
    tempIndx = find_amisr_exp(superMag.time(iStorm),amisr.startTime, amisr.endTime);
    numExp(iStorm) = length(tempIndx);
    amisrIndx(iStorm)=tempIndx(1);
end
superMag.expID = repmat(string("nan"),1,length(superMag.stormID));
superMag.expName = repmat(string("nan"),1,length(superMag.stormID));
superMag.status = repmat(string("nan"),1,length(superMag.stormID));
superMag.startTime = nan(1,length(superMag.stormID));
superMag.endTime = nan(1,length(superMag.stormID));
superMag.expBC = false(1,length(superMag.stormID));
superMag.numberOfSimultaneousExp = nan(1,length(superMag.stormID));

superMag.expID(~isnan(amisrIndx))=amisr.expId(amisrIndx(~isnan(amisrIndx)));
superMag.expName(~isnan(amisrIndx))=amisr.expName(amisrIndx(~isnan(amisrIndx)));
superMag.status(~isnan(amisrIndx))=amisr.status(amisrIndx(~isnan(amisrIndx)));
superMag.startTime(~isnan(amisrIndx))=amisr.startTime(amisrIndx(~isnan(amisrIndx)));
superMag.endTime(~isnan(amisrIndx))=amisr.endTime(amisrIndx(~isnan(amisrIndx)));
superMag.expBC(~isnan(amisrIndx))=bcFilterIndx(amisrIndx(~isnan(amisrIndx)));
superMag.numberOfSimultaneousExp(~isnan(amisrIndx))=numExp(amisrIndx(~isnan(amisrIndx)));


%% Create a table
T = table(superMag.datetime,superMag.stormID,...
    superMag.AE, superMag.mlat,...
    superMag.mlt, superMag.pfisrMlt,...
    superMag.expID',superMag.expName',...
    datetime(superMag.startTime,'ConvertFrom','datenum')',...
    datetime(superMag.endTime,'ConvertFrom','datenum')',...
    superMag.status',...
    superMag.expBC',...
    'VariableNames',{'Time','stormID','AE','MLAT','MLT','PFISR_MLT',...
    'PFISR_ExpID','PFISR_ExpName','PFISR_startTime','PFISR_endTime',...
    'PFISR_ExpStatus','BarkerCode'});
%%
% disp('Table of Barker Code Experiments during a SuperMag substorm');
% T(T.BarkerCode,:)
%%
% disp('Table of PFISR Experiments during a SuperMag substorm');
% T(~strcmp(T.PFISR_ExpID,"nan"),:)

%% Finding if DASC is ON during a particular substorm with PFISR ON
substormIndx = 1:1:length(superMag.time);
cIndx = substormIndx(closestSubstormIndx);
Tdasc = read_h5_data(dascFileStr);
[timeStamp, wavelength] = restructure_DASC_table_to_time_array(Tdasc);

timeStamp = datenum(timeStamp);
superMagCloseStorm.time = datetime(superMag.time(cIndx),'ConvertFrom','datenum');
tempStart = datenum(dateshift(superMagCloseStorm.time,'start','day'));
tempEnd = datenum(dateshift(superMagCloseStorm.time,'end','day'));
superMagCloseStorm.DASC_timeMin = repmat(NaT,length(cIndx),1);
superMagCloseStorm.DASC_timeMax = repmat(NaT,length(cIndx),1);
superMagCloseStorm.DASC_wavelength = repmat({"nan"},length(cIndx),1);

Ft = griddedInterpolant(timeStamp,1:numel(timeStamp),'nearest','nearest');
superMagCloseStorm.DASC_timeMin = datetime(timeStamp(Ft(tempStart)),'ConvertFrom','datenum');
superMagCloseStorm.DASC_timeMax = datetime(timeStamp(Ft(tempEnd)),'ConvertFrom','datenum');
tmIndex = superMagCloseStorm.DASC_timeMin==superMagCloseStorm.DASC_timeMax;
superMagCloseStorm.DASC_timeMin(tmIndex) = NaT;
superMagCloseStorm.DASC_timeMax(tmIndex) = NaT;

for cStorm = 1:1:length(cIndx)
    if Ft(tempEnd(cStorm)) ~= Ft(tempStart(cStorm))
    superMagCloseStorm.DASC_wavelength(cStorm,1) = {unique(wavelength(Ft(tempStart(cStorm)):Ft(tempEnd(cStorm))))};
    end
end

 %% Substorms that are close to PFISR
 T1 = table(superMag.datetime(closestSubstormIndx),superMag.stormID(closestSubstormIndx),...
    superMag.AE(closestSubstormIndx), superMag.mlat(closestSubstormIndx),...
    superMag.mlt(closestSubstormIndx), superMag.pfisrMlt(closestSubstormIndx),...
    superMag.expID(closestSubstormIndx)',superMag.expName(closestSubstormIndx)',...
    datetime(superMag.startTime(closestSubstormIndx)','ConvertFrom','datenum'),...
    datetime(superMag.endTime(closestSubstormIndx)','ConvertFrom','datenum'),...
    superMag.status(closestSubstormIndx)',...
    superMag.expBC(closestSubstormIndx)',...
    superMagCloseStorm.DASC_timeMin,...
    superMagCloseStorm.DASC_timeMax,...
    superMagCloseStorm.DASC_wavelength,...
    'VariableNames',...
    {'Time','stormID','AE','MLAT','MLT','PFISR_MLT','PFISR_ExpID',...
    'PFISR_ExpName','PFISR_startTime','PFISR_endTime',...
    'PFISR_ExpStatus','BarkerCode',...
    'DASC_TimeMin','DASC_TimeMax','DASC_Wavelength'});

end

function [time, wavelength, url, wavelengthStr] = get_DASC_times_during_substorm(...
    dascTimeStamps, wavelengths, substormTime, growthDuration, expansionDuration)
    
    %substormTime - datetime
    if nargin<5
        expansionDuration = 1.0; %hr
    end
    
    if nargin<4
        growthDuration = 2.0; %hr
    end
    
    tempStart = substormTime - hours(growthDuration);
    tempEnd = substormTime + hours(expansionDuration);
    indxtStamp = dascTimeStamps<=tempEnd & dascTimeStamps>=tempStart;
    time = dascTimeStamps(indxtStamp);
    wavelength = wavelengths(indxtStamp);
    wavelengthStr = num2str(wavelength,'%04.f');
    url = create_DASC_url(time,wavelength);
end

function data=load_ascii_files(loadFile, format, headerlines)
if nargin<3
    headerlines = 1;
end
fileID = fopen(loadFile,'r');
data = textscan(fileID, format, 'headerlines', headerlines);
fclose(fileID);
end

function superMag = extract_superMag_data(loadFile)
format ='%4f %2f %2f %2f %2f %5.2f %5.2f ';
tempData = load_ascii_files(loadFile, format, 69);
superMag.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
superMag.time = datenum(datestr(superMag.datetime));
superMag.mlat = tempData{6};
superMag.mlt = tempData{7};
end

function amisr = extract_amisr_data(loadFile)
    table = read_h5_data(loadFile);
    amisr.expId = string(table.Data{2});
    amisr.expName = string(table.Data{3});
    amisr.status = string(table.Data{4});
    amisr.startTime = unixtime2matlab(table.Data{5});
    amisr.endTime = unixtime2matlab(table.Data{1});
end

function [Lm,MLT] = get_pfisr_magnetic_coordinates(time,maginput,GDZ,magFieldNo)
    [Lm,~,~,~,~,MLT] = onera_desp_lib_make_lstar(magFieldNo,[0,0,0,0,0],0,time,GDZ(3),GDZ(1),GDZ(2),maginput);
%     tup=py.aacgmv2.wrapper.get_aacgm_coord(GDZ(1), GDZ(2), GDZ(3), time, 'TRACE');
%     MLAT = double(py.array.array('d',py.numpy.nditer(tup{1})));
%     MLON = double(py.array.array('d',py.numpy.nditer(tup{2})));
%     MLT_AACGM = double(py.array.array('d',py.numpy.nditer(tup{3})));
    Lm = abs(Lm);
end

function [amisrIndx] = find_amisr_exp(time, startTimeArr, endTimeArr)
% Finds the amisr experiment indx
    amisrIndx = find(time>startTimeArr & time<endTimeArr);
    if isempty(amisrIndx)
        amisrIndx=nan;
    end
    
end

function expArr = barker_coded_experiments()
expArr =["GenPINOT_PulsatingAurora_TN30          ";
    "Inspire_v01                            ";
    "Kelley01                               ";
    "MSWinds23                              ";
    "MSWinds23_dt013                        ";
    "MSWinds23_3dt                            ";
    "MSWinds23m                               ";
    "MSWinds23hr                              ";
    "MSWinds21                                ";
    "MSWinds26.v03                          ";
    "Semeter01                              ";
    "Sporadic01                             ";
    "Sporadic02                             ";
    "Sporadic03                             ";
    "Sporadic04                             ";
    "Sporadic14                             ";
    "Sporadic15                             ";
    "Sporadic15_3dt                         ";
    "ThemisD1.v01                             "
    "MSWinds27.v01                            "];
end
 
function [timeStamp, wavelength] = restructure_DASC_table_to_time_array(T)
    wavelength = string();
    timeStamp=cell2mat(T.Data(strcmp(string(T.Name),"time")));
    waveCode=T.Data(strcmp(string(T.Name),"wavelengthCode"));
    waveLength=T.Data(strcmp(string(T.Name),"wavelength"));
    for i = 0:1:length(waveCode)-1
        if i==0
            wavelength(end:end+numel(waveLength{i+1})-1,1) = string(waveCode{i+1}(waveLength{i+1},1));
        else
            wavelength(end+1:end+numel(waveLength{i+1}),1) = string(waveCode{i+1}(waveLength{i+1},1));
        end
    end
    [timeStamp,I] = sort(timeStamp);
    [timeStamp, w] = unique(timeStamp, 'stable');
    wavelength = wavelength(I);
    wavelength = wavelength(w);
    timeStamp = datetime(unix_to_matlab_time(timeStamp),'ConvertFrom','datenum');
end

function dMLT = absDiffMLT(a,b)
    % Source: https://stackoverflow.com/questions/32276369/calculating-absolute-differences-between-two-angles 
    normMLT = mod(a-b,24);
    dMLT = min(24-normMLT, normMLT);
end
