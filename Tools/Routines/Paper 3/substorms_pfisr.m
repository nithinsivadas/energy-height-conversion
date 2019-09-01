%% Substorms in the vicinity of PFISR
% Find SuperMag Substorms and PFISR conjunctions

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

load([storeDir,'table_of_substorms_as_input.mat']);

% Temporary file

rawFileListStr = select_PFISR_raw_data(storeDir);

%% Loading database
timeMinStr = "01 Dec 2006";
timeMaxStr = "31 Jul 2019";

% Table of substorms where DASC is ON
T2 = T1(~strcmp(T1.DASC_Wavelength,'nan') & T1.Time<T1.DASC_TimeMax & T1.Time>T1.DASC_TimeMin,:);
% Table of substorms with an additional constraint of an operating PFISR barker code mode
T3 = T2(T2.BarkerCode,:);
% Table of substorms within an additional constraing of being within the user defined dates
T4 = T3(T3.Time>=datetime(datestr(timeMinStr)) & T3.Time<=datetime(datestr(timeMaxStr)),:);

%%
for j=1:1:length(rawFileListStr)
    tempStr = strsplit(rawFileListStr(j),filesep);
    tempStr = strsplit(tempStr(end),'_');
    expID(j,:) = tempStr(1);
end

%%
for i = 1:1:length(rawFileListStr)
    stormIDIndx = find(strcmp(deblank(T4.PFISR_ExpID),expID(i)));
    
    for k=1:1:length(stormIDIndx)
    % Can be a batch file if we need to
    stormTime = T4.Time(stormIDIndx(k));
    minTime = stormTime-hours(2);
    maxTime = stormTime+hours(1);
    outputPrefix = num2str(T4.stormID(stormIDIndx(k)));
    write_ne_to_h5_v2(outputPrefix, outputh5Suffix, rawFileListStr(i),storeDir,...
        minTime, maxTime, stormTime);
    end
end

%% Conductivity and Energy
fileNameList = struct2cell(dir([storeDir,'*_pfisrData.h5']));
filePathStr = strcat(storeDir,string(fileNameList(1,:)'));

energyBin = logspace(3,6,25)';

for i=1:1:length(filePathStr)
    fileName = filePathStr(i);
    [data1] = calculate_conductivity(fileName);
    write_conductivity(data1,fileName);
    [data2] = calculate_energy(fileName,energyBin);
    write_energy(data2,fileName);
end

%% Functions
function [data] = calculate_energy(fileName,energyBin)
   
   stepsHour = hours(1); %Calculates conductivity and writes it in steps of x hours
   row = @(x,y) strcmp(string(x.Path),y);
   in = read_h5_data(fileName);
   alt = in.Data{row(in,'/alt')}';
   electronDensity = in.Data{row(in,'/inputData/Ne')};
   dNeFrac = in.Data{row(in,'/inputData/dNeFrac')}; % dNeFrac is actually dNe
   latitude = in.Data{row(in,'/site/latitude')};
   longitude = in.Data{row(in,'/site/longitude')};
   time = in.Data{row(in,'/time')};
   ttime = datetime(time,'ConvertFrom','datenum');

    
%    mode = 0;
%    outputs = [];
%    setPlotConductivity = false;
   
   hourArr = ttime(1):stepsHour:ttime(end);
   nhour = length(hourArr);

   for iHour = 1:1:nhour
       minTimeIndx = find_time(time,datenum(hourArr(iHour)));
       if iHour ~=nhour 
        maxTimeIndx = find_time(time,datenum(hourArr(iHour+1)));
       else
        maxTimeIndx = length(time);
       end
       medTimeIndx = round((minTimeIndx + maxTimeIndx)/2);
       lat1=latitude*ones(size(alt));
       lon1=longitude*ones(size(alt));
       A = get_energy_dep_matrix(alt,energyBin,lat1,lon1,time(medTimeIndx));
       
       timeIndx = minTimeIndx:1:maxTimeIndx;
       
       if length(timeIndx)==1
       data.energyFlux(timeIndx,:) = nan;
        data.denergyFlux(timeIndx,:) = nan;
        data.MSE(timeIndx,:) = nan;
        data.qInverted(timeIndx,:) = nan;
        data.q(timeIndx,:) = nan;
        data.dq(timeIndx,:) = nan;
       else
       [dq, q] = get_error_in_q(electronDensity(timeIndx,:)',...
            dNeFrac(timeIndx,:)',alt,time(timeIndx)',2);
        dataInv = get_inverted_flux(q,dq,time(timeIndx),alt,energyBin,A);
        data.energyFlux(timeIndx,:) = dataInv.energyFlux';
        data.denergyFlux(timeIndx,:) = dataInv.dEnergyFlux';
        data.MSE(timeIndx,:) = dataInv.MSE';
        data.qInverted(timeIndx,:) = dataInv.qInverted;
        data.q(timeIndx,:) = q;
        data.dq(timeIndx,:) = dq;       
       end
        
        
        data.time(timeIndx) = time(timeIndx);
       
   end
   data.energyBin = energyBin;
   data.altitude = alt;
end

function write_energy(data,fileName)
    write_h5_dataset(fileName,'/energy/energyFlux',data.energyFlux,1);
    write_h5_dataset(fileName,'/energy/dEnergyFlux',data.denergyFlux,1);
    write_h5_dataset(fileName,'/energy/energyBin',data.energyBin,0);
    write_h5_dataset(fileName,'/energy/MSE',data.MSE,1);
    write_h5_dataset(fileName,'/energy/q',data.q,1);
    write_h5_dataset(fileName,'/energy/dq',data.dq,1);
    write_h5_dataset(fileName,'/energy/qInverted',data.qInverted,1);
end


function [data] = calculate_conductivity(fileName)
   
   stepsHour = hours(1); %Calculates conductivity and writes it in steps of x hours
   row = @(x,y) strcmp(string(x.Path),y);
   in = read_h5_data(fileName);
   alt = in.Data{row(in,'/alt')}';
   electronDensity = in.Data{row(in,'/inputData/Ne')};
   latitude = in.Data{row(in,'/site/latitude')};
   longitude = in.Data{row(in,'/site/longitude')};
   time = in.Data{row(in,'/time')};
   ttime = datetime(time,'ConvertFrom','datenum');

    
   mode = 0;
   outputs = [];
   setPlotConductivity = false;
   
   hourArr = ttime(1):stepsHour:ttime(end);
   nhour = length(hourArr);
   alt1=alt(:)';
   lat1=latitude*ones(size(alt1));
   lon1=longitude*ones(size(alt1));
   for iHour = 1:1:nhour
       minTimeIndx = find_time(time,datenum(hourArr(iHour)));
       if iHour ~=nhour 
        maxTimeIndx = find_time(time,datenum(hourArr(iHour+1)));
       else
        maxTimeIndx = length(time);
       end
       medTimeIndx = round((minTimeIndx + maxTimeIndx)/2);
       iriData = iri2016f90(time(medTimeIndx), alt, latitude, longitude);
       time1=time(medTimeIndx)*ones(size(alt1));
       msisData = msis_irbem(time1, [alt1',lat1',lon1']);
       for iTime = minTimeIndx:1:maxTimeIndx
        [tempData] = get_conductivity_v2( alt, electronDensity(iTime,:)',...
        latitude, longitude, time(iTime), mode, outputs,setPlotConductivity,...
        iriData, msisData);
        data.pedersonConductivity(iTime,:) = tempData.pedersonConductivity;
        data.hallConductivity(iTime,:) = tempData.hallConductivity;
        data.time(iTime) = time(iTime);
       end
   end
   data.altitude = tempData.altitude;
end

function write_conductivity(data,fileName)
    write_h5_dataset(fileName,'/conductivity/hall',real(data.hallConductivity),1);
    write_h5_dataset(fileName,'/conductivity/pederson',real(data.pedersonConductivity),1);
    write_h5_dataset(fileName,'/conductivity/ihall',imag(data.hallConductivity),1);
    write_h5_dataset(fileName,'/conductivity/ipederson',imag(data.pedersonConductivity),1);
end


function write_ne_to_h5_v2(outputPrefix,outputh5Suffix,rawFileStr,storeDir,...
    minTime,maxTime,stormTime)
 
    % specific data structure written to h5 file
    data = extract_ne_from_raw_h5(rawFileStr);
    minTimeIndx = find_time(data.time,datestr(minTime));
    maxTimeIndx = find_time(data.time,datestr(maxTime));
    tIndx = minTimeIndx:1:maxTimeIndx;
    tempStr = strsplit(rawFileStr,filesep);
    tempStr = strsplit(tempStr(end),'_');
    expID = tempStr(1);
    outputh5Str = strcat(outputPrefix,'_',datestr(stormTime,'YYYYmmdd'),'_',outputh5Suffix);
    outputh5Str = strcat(storeDir,outputh5Str);
    
    % Writing Electron Density
    write_h5_dataset(outputh5Str,'/inputData/Ne',data.electronDensity(tIndx,:),1,true);
    try
    write_h5_dataset(outputh5Str,'/inputData/dNeFrac',data.dNeFrac(tIndx,:),1,true);
    catch ME
        warning('No dNeFrac to write.');
    end
    
    % Writing Coordinates
    write_h5_dataset(outputh5Str,'/time',data.time(tIndx,:),1,true);
    write_h5_dataset(outputh5Str,'/lat',data.lat,0,true);
    write_h5_dataset(outputh5Str,'/lon',data.lon,0,true);
    write_h5_dataset(outputh5Str,'/az',data.az,0,true);
    write_h5_dataset(outputh5Str,'/el',data.el,0,true);
    write_h5_dataset(outputh5Str,'/range',data.range,0,true);
    write_h5_dataset(outputh5Str,'/alt',data.alt,0,true);
    
    write_h5_dataset(outputh5Str,'/site/latitude',data.site.latitude,0,true);
    write_h5_dataset(outputh5Str,'/site/longitude',data.site.longitude,0,true);
    write_h5_dataset(outputh5Str,'/site/altitude',data.site.altitude,0,true);
    h5_create_writestr(outputh5Str,'/inputFile',rawFileStr);
    h5_create_writestr(outputh5Str,'/expID',expID);
%     write_h5_dataset(outputh5Str,'/inputFile',rawFileStr,0,true);
    
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
    write_h5_dataset(outputh5Str,'/alt',data.alt,0,true);
    
    write_h5_dataset(outputh5Str,'/site/latitude',data.site.latitude,0,true);
    write_h5_dataset(outputh5Str,'/site/longitude',data.site.longitude,0,true);
    write_h5_dataset(outputh5Str,'/site/altitude',data.site.altitude,0,true);
    
end

function data = extract_ne_from_raw_h5(rawFileStr)
    
    tempData = read_amisr(rawFileStr);
    data.electronDensity = permute(tempData.electronDensity(:,tempData.magBeamNo,:),[3 1 2]);
    data.altitude = tempData.altitude(:,tempData.magBeamNo)';
    data.time = tempData.time(1,:)';
    data.dNeFrac = permute(tempData.dNeFrac(:,tempData.magBeamNo,:),[3 1 2]);
    data.range = tempData.range(:,tempData.magBeamNo)';
    data.az = tempData.az(:,tempData.magBeamNo)';
    data.el = tempData.el(:,tempData.magBeamNo)';
    [data.lat,data.lon,data.alt] = aer2geodetic(data.az,data.el,data.range,...
        tempData.site.latitude,tempData.site.longitude,tempData.site.altitude/1000,wgs84Ellipsoid('km'));
    data.site = tempData.site;
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



