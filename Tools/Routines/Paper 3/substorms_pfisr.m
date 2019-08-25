%% Substorms in the vicinity of PFISR
% Find SuperMag Substorms and PFISR conjunctions
clear all;

%% Initialization
if strcmp(get_computer_name,'nithin-carbon')
    dataDir = 'G:\My Drive\Research\Projects\Data\';
    storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\';
elseif strcmp(get_computer_name,'aurora1-optiplex-780')
    dataDir = '/media/nithin/PFISR_002_006/Nithin/Data/';
    storeDir = '/media/nithin/PFISR_002_006/Nithin/Data/Paper_3/';
elseif strcmp(get_computer_name,'scc-lite')
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/scratch/nithin/Data/Paper_3/';
else
%     error(['Not configured for this computer: ',get_computer_name]);
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/projectnb/semetergrp/nithin/Data/Paper_3/';
end

outputAMISRFileStr = 'amisrWebDatabase.h5';
amisrDatabaseStr = [dataDir,outputAMISRFileStr];
superMagFileStr = [storeDir,'substorms_superMag_20190530.txt'];
dascFileStr = [storeDir,'dascDatabase.h5'];
omniFileStr = [dataDir,'omni.h5'];


outputh5Suffix = 'pfisrData.h5';

% Substorms at PFISR [IMPORTANT]
Dmlt = 2;
Dmlat = 1;

% Poker Flat geolocation
pkrGLAT = 65.126;
pkrGLON = -147.47;
pkrh0=0.693;

% tic 
% [T,T1, timeStamp, wavelength] = substorm_create_table(dataDir,outputAMISRFileStr,amisrDatabaseStr,...
%   superMagFileStr, dascFileStr, omniFileStr,...
%   timeMinStr, timeMaxStr, Dmlt, Dmlat, pkrGLAT, pkrGLON, pkrh0);
% toc 

% load([storeDir,'table_of_substorms_as_input.mat']);

% Temporary file

rawFileListStr = select_PFISR_raw_data(storeDir);

for i = 1:1:length(rawFileListStr)
    % Can be a batch file if we need to
    write_ne_to_h5(outputh5Suffix, rawFileListStr(i),storeDir);
end

function write_ne_to_h5(outputh5Suffix,rawFileStr,storeDir)
 
    % specific data structure written to h5 file
    data = extract_ne_from_raw_h5(rawFileStr);
    tempStr = strsplit(rawFileStr,filesep);
    tempStr = strsplit(tempStr(end),'_');
    expID = tempStr(1);
    outputh5Str = strcat(expID,'_',outputh5Suffix);
    outputh5Str = strcat(storeDir,outputh5Str);
    
    % Writing Electron Density
    write_h5_dataset(outputh5Str,'/inputData/Ne',data.electronDensity,1,true);
    try
    write_h5_dataset(outputh5Str,'/inputData/dNeFrac',data.dNeFrac,1,true);
    catch ME
        warning('No dNeFrac to write.');
    end
    
    % Writing Coordinates
    write_h5_dataset(outputh5Str,'/time',data.time,1,true);
    write_h5_dataset(outputh5Str,'/lat',data.lat,0,true);
    write_h5_dataset(outputh5Str,'/lon',data.lon,0,true);
    write_h5_dataset(outputh5Str,'/az',data.az,0,true);
    write_h5_dataset(outputh5Str,'/el',data.el,0,true);
    write_h5_dataset(outputh5Str,'/range',data.range,0,true);
    
end

function data = extract_ne_from_raw_h5(rawFileStr)
    
    tempData = read_amisr(rawFileStr);
    data.electronDensity = permute(tempData.electronDensity(:,tempData.magBeamNo,:),[3 1 2]);
    data.altitude = tempData.altitude(:,tempData.magBeamNo)';
    data.time = permute(tempData.electronDensity(:,tempData.magBeamNo,:),[3 1 2]);
    data.dNeFrac = permute(tempData.dNeFrac(:,tempData.magBeamNo,:),[3 1 2]);
    data.range = tempData.range(:,tempData.magBeamNo)';
    data.az = tempData.az(:,tempData.magBeamNo)';
    data.el = tempData.el(:,tempData.magBeamNo)';
    [data.lat,data.lon,~] = aer2geodetic(data.az,data.el,data.range,...
        tempData.site.latitude,tempData.site.longitude,tempData.site.altitude/1000,wgs84Ellipsoid('km'));
end

function rawFileListStr = select_PFISR_raw_data(storeDir)
        
    dataDir = strcat(storeDir,'PFISR',filesep);
   
    if isfolder(dataDir)
        rawFileList = struct2cell(dir([dataDir,'*.h5']))'; 
    else
        error(strcat(dataDir, ' does not exist'));
    end
    
    rawFileListStr = string(strcat(rawFileList(:,2),filesep,rawFileList(:,1)));
end



