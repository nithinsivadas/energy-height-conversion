%% PFISR Processing based on storm list
% Find substorms and pfisr conjunctions
% Calculate and store data - electron density, conductivity and energy flux
% from already existing pfisr data in folder

if strcmp(get_computer_name,'nithin-carbon')
    dataDir = 'G:\My Drive\Research\Projects\Data\';
    storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\Processed\';
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

substormListFileStr = [storeDir,'sophie_output_90.asc'];
substormListFileType = 2; % 1 -> superMag database, 2 -> Colin Forsyth's data
omniFile = 'G:\My Drive\Research\Projects\Data\omni.h5';
workDir = 'G:\My Drive\Research\Projects\Paper 3\Data\WorkDir';

outputh5Suffix = 'pfisrData.h5';
minAlt = 40;

timeMinStr = "01 Dec 2006";
timeMaxStr = "31 Dec 2019";
[T, superMag] = substorm_create_table(substormListFileStr, substormListFileType,...
  timeMinStr, timeMaxStr);

%%
rawFileListStr = select_PFISR_raw_data(storeDir);

T3 = T(T.BarkerCode,:);

T4 = T3(T3.Time>=datetime(datestr(timeMinStr)) & T3.Time<=datetime(datestr(timeMaxStr)),:);

%%
for j=1:1:length(rawFileListStr)
    tempStr = strsplit(rawFileListStr(j),filesep);
    tempStr = strsplit(tempStr(end),'_');
    expID(j,:) = tempStr(1);
end

%% Generate processed PFISR files
dn = 1./length(rawFileListStr);
multiWaitbar('Extract AMISR',0);

    for i=1:1:1
%     for i=1:1:length(rawFileListStr)
        % Check if the raw pfisr data has a substorm during it
        rawFileQualify(i)=sum(strcmp(deblank(T4.PFISR_ExpID),expID(i)))>=1; 
        
        if rawFileQualify(i)
            outputPrefix = expID(i);
            outputh5Suffix = 'pfisrData.h5';
            write_ne_to_h5_v3(outputPrefix,outputh5Suffix, rawFileListStr(i),storeDir,...
        minAlt);
            
        end
        
        multiWaitbar('Extract AMISR','Increment',dn);    
    end

%% Conductivity and Energy
fileNameList = struct2cell(dir([storeDir,'*_pfisrData.h5']));
filePathStr = strcat(storeDir,string(fileNameList(1,:)'));

energyBin = logspace(3,6,25)';
multiWaitbar('Calculate Conductivity',0);
multiWaitbar('Calculate Energy',0);
dn = 1./length(filePathStr);
for i=1:1:1
% for i=64:1:length(filePathStr)
    fileName = filePathStr(i);
    [data1] = calculate_conductivity(fileName);
    write_conductivity(data1,fileName);
    multiWaitbar('Calculate Conductivity','Increment',dn);
    [data2] = calculate_energy(fileName,energyBin);
    write_energy(data2,fileName);
    multiWaitbar('Calculate Energy','Increment',dn);
    
end
multiWaitbar('Close All');
disp('Complete');

%%

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

function write_ne_to_h5_v3(outputPrefix,outputh5Suffix,rawFileStr,storeDir,...
    minAlt)
 
    % specific data structure written to h5 file
    data = extract_ne_from_raw_h5(rawFileStr,minAlt);
    minTimeIndx = find_time(data.time,data.time(1));
    maxTimeIndx = find_time(data.time,data.time(end));
    tIndx = minTimeIndx:1:maxTimeIndx;
    tempStr = strsplit(rawFileStr,filesep);
    tempStr = strsplit(tempStr(end),'_');
    expID = tempStr(1);
    outputh5Str = strcat(outputPrefix,'_',outputh5Suffix);
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

function data = extract_ne_from_raw_h5(rawFileStr, minAlt)
    
    tempData = read_amisr(rawFileStr);
    data.electronDensity = permute(tempData.electronDensity(:,tempData.magBeamNo,:),[3 1 2]);
    data.altitude = tempData.altitude(:,tempData.magBeamNo)';
    data.time = tempData.time(1,:)';
    try
    data.dNeFrac = permute(tempData.dNeFrac(:,tempData.magBeamNo,:),[3 1 2]);
    catch ME
    end
    data.range = tempData.range(:,tempData.magBeamNo)';
    data.az = tempData.az(:,tempData.magBeamNo)';
    data.el = tempData.el(:,tempData.magBeamNo)';
    [data.lat,data.lon,data.alt] = aer2geodetic(data.az,data.el,data.range,...
        tempData.site.latitude,tempData.site.longitude,tempData.site.altitude/1000,wgs84Ellipsoid('km'));
    data.site = tempData.site;
    
        %% Trim lower limit of altitude appropriately
    altIndx = find_altitude(data.alt,minAlt); % km
    indx = altIndx:1:length(data.alt);
    
    data.electronDensity=data.electronDensity(:,indx);
    data.altitude=data.altitude(1,indx);
    data.dNeFrac=data.dNeFrac(:,indx);
    data.range = data.range(1,indx);
    data.az = data.az(1,indx);
    data.el = data.el(1,indx);
    data.lat = data.lat(1,indx);
    data.lon = data.lon(1,indx);
    data.alt = data.alt(1,indx);
    
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


function amisr = extract_amisr_data(loadFile)
    table = read_h5_data(loadFile);
    amisr.expId = string(table.Data{2});
    amisr.expName = string(table.Data{3});
    amisr.status = string(table.Data{4});
    amisr.startTime = unixtime2matlab(table.Data{5});
    amisr.endTime = unixtime2matlab(table.Data{1});
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

% A function that extracts data from substorm/storm/SCM lists
function T1 = extract_data(loadFile, ftype)
    
    if ftype == 1
        format ='%4f %2f %2f %2f %2f %5.2f %5.2f ';
        tempData = load_ascii_files(loadFile, format, 69);
        superMag.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        superMag.time = datenum(datestr(superMag.datetime));
        superMag.mlat = tempData{6};
        superMag.mlt = tempData{7};
    end

    if ftype == 2 %Forsyth SOPHIE substorm phases data
        format = '%4f/%2f/%2f-%2f:%2f:%2f %u %u %5.2f %5.2f';
        tempData = load_ascii_files(loadFile, format, 16);

        T.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},tempData{6});
        T.time = datenum(datestr(T.datetime));
        T.phase = tempData{7};
        T.flag = tempData{8};
        T.mlt = tempData{9};
        T.mlat = tempData{10};
    end
    
    if ftype == 3 %Walach Storm data
        format = '%4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f';
        tempData = load_ascii_files(loadFile, format, 6);

        T.datetimeI = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        T.datetimeM = datetime(tempData{1+5},tempData{2+5},tempData{3+5},tempData{4+5},tempData{5+5},zeros(size(tempData{5+5})));
        T.datetimeR = datetime(tempData{1+5+5},tempData{2+5+5},tempData{3+5+5},tempData{4+5+5},tempData{5+5+5},zeros(size(tempData{5+5+5})));
        T.datetimeE = datetime(tempData{1+5+5+5},tempData{2+5+5+5},tempData{3+5+5+5},tempData{4+5+5+5},tempData{5+5+5+5},zeros(size(tempData{5+5+5+5})));
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeM = datenum(datestr(T.datetimeM));
        T.timeR = datenum(datestr(T.datetimeR));
        T.timeE = datenum(datestr(T.datetimeE));
        T.symH_min = tempData{21};
    end
    
     if ftype == 4 %Walach SCM data with preceeding pre without preceeding substorms
        format = '%4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f';
        tempData = load_ascii_files(loadFile, format, 15);
        T.datetimeI = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        T.datetimeE = datetime(tempData{1+5},tempData{2+5},tempData{3+5},tempData{4+5},tempData{5+5},zeros(size(tempData{5+5})));
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeE = datenum(datestr(T.datetimeE));
        
     end
    
     if ftype == 5 %Luisa, current sheet scattering
        format = ['%s %4f-%2f-%2f/%2f:%2f:%2f.%3f \t %4f-%2f-%2f/%2f:%2f:%2f.%3f',repmat(' %*s',1,44)];
        tempData = load_ascii_files(loadFile, format, 0);
        T.sc = tempData{1};
        T.datetimeI = datetime(tempData{2},tempData{3},tempData{4},tempData{5},tempData{6},tempData{7},tempData{8});
        T.datetimeE = datetime(tempData{1+8},tempData{2+8},tempData{3+8},tempData{4+8},tempData{5+8},tempData{6+8},tempData{7+8});
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeE = datenum(datestr(T.datetimeE));
        
    end
    
    T1 = struct2table(T);
end

function [T, superMag] = substorm_create_table(superMagFileStr, superMagFileType,...
  timeMinStr, timeMaxStr)

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
% dascFileStr = [storeDir,'dascDatabase.h5'];
omniFileStr = [dataDir,'omni.h5'];

% Substorms at PFISR [IMPORTANT] : Range of closeness to PFISR
% Dmlt = 2;
% Dmlat = 50; %desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat

% Poker Flat geolocation
pkrGLAT = 65.126;
pkrGLON = -147.47;
pkrh0=0.693;

% Loading database

% Create amisr database
if ~isfile(amisrDatabaseStr)
    status = create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,[dataDir,outputAMISRFileStr]);
end
%
% Load supermag database
superMag1 = extract_data(superMagFileStr,superMagFileType);
superMag = superMag1(superMag1.phase==2,:);
superMag.stormID = (1:length(superMag.time))';

for i=1:1:length(superMag.stormID)
    if i==1
        superMag.growthPhaseStart(i) =nan;
    else
        try 
        superMag.growthPhaseStart(i)=superMag1{superMag1.time > superMag.time(i-1) & superMag1.time < superMag.time(i) & superMag1.phase==1,2};
        catch
            superMag.growthPhaseStart(i)=nan;
        end
    end
    
    if i==length(superMag.stormID)
        superMag.recoveryPhaseStart(i) = nan;
        superMag.recoveryPhaseEnd(i) = nan;
    else
        try
        superMag.recoveryPhaseStart(i)=superMag1{superMag1.time > superMag.time(i) & superMag1.time < superMag.time(i+1) & superMag1.phase==3,2};
        catch
            superMag.recoveryPhaseStart(i) = nan;
        end
        
        try
            superMag.recoveryPhaseEnd(i)=superMag1{superMag1.time > superMag.time(i) & superMag1.time < superMag.time(i+1) & superMag1.phase~=3 ,2};
        catch
            superMag.recoveryPhaseEnd(i)=nan;
        end
        
    end
        
end

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
% deltaMLT = absDiffMLT(superMag.pfisrMlt,superMag.mlt);
% desiredMLTIndx = abs(deltaMLT)<Dmlt;
% desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat;
% closestSubstormIndx = desiredMLTIndx & desiredMLATIndx; 

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
superMag.expID = repmat(string("nan"),1,length(superMag.stormID))';
superMag.expName = repmat(string("nan"),1,length(superMag.stormID))';
superMag.status = repmat(string("nan"),1,length(superMag.stormID))';
superMag.startTime = nan(1,length(superMag.stormID))';
superMag.endTime = nan(1,length(superMag.stormID))';
superMag.expBC = false(1,length(superMag.stormID))';
superMag.numberOfSimultaneousExp = nan(1,length(superMag.stormID))';

superMag.expID(~isnan(amisrIndx))=amisr.expId(amisrIndx(~isnan(amisrIndx)));
superMag.expName(~isnan(amisrIndx))=amisr.expName(amisrIndx(~isnan(amisrIndx)));
superMag.status(~isnan(amisrIndx))=amisr.status(amisrIndx(~isnan(amisrIndx)));
superMag.startTime(~isnan(amisrIndx))=amisr.startTime(amisrIndx(~isnan(amisrIndx)));
superMag.endTime(~isnan(amisrIndx))=amisr.endTime(amisrIndx(~isnan(amisrIndx)));
superMag.expBC(~isnan(amisrIndx))=bcFilterIndx(amisrIndx(~isnan(amisrIndx)));
superMag.numberOfSimultaneousExp(~isnan(amisrIndx))=numExp(amisrIndx(~isnan(amisrIndx)));


%% Create a table
T = table(superMag.datetime,superMag.stormID,...
    datetime(superMag.growthPhaseStart,'ConvertFrom','datenum'),...
    datetime(superMag.recoveryPhaseStart,'ConvertFrom','datenum'),...
    datetime(superMag.recoveryPhaseEnd,'ConvertFrom','datenum'),...
    superMag.AE, superMag.mlat,...
    superMag.mlt, superMag.pfisrMlt,...
    superMag.expID,superMag.expName,...
    datetime(superMag.startTime,'ConvertFrom','datenum'),...
    datetime(superMag.endTime,'ConvertFrom','datenum'),...
    superMag.status,...
    superMag.expBC,...
    'VariableNames',{'Time','stormID','GPStart','RPStart','RPEnd','AE','MLAT','MLT','PFISR_MLT',...
    'PFISR_ExpID','PFISR_ExpName','PFISR_startTime','PFISR_endTime',...
    'PFISR_ExpStatus','BarkerCode'});

end