function [omniData, matFilePath] = process_omni_data(dateStr,dataStoreDir)
%% process_omni_cdf.m Processes hourly and minute-wise omni data and stores 
%                     it as .mat files, under each event date folder. 
%                     Can be used as input for magnetic field models
%                     through IRBEM. 
%--------------------------------------------------------------------------
% Input
%------
% dataStr - [dd mmm yyyy] A string indicating date for which you would like
%           to download the omni data
%           Default:'26 Mar 2008'
% storeDir- A string indicating the data storage directory in '...\github\'
%           It will store the hourly and 1 min data in ~yyyymmdd\omni\1min\
%           Default: '~\github\LargeFiles\'
%--------------------------------------------------------------------------
% Output
%-------
% Stores the .mat files in dataStoreDir.
% omniData -> hourly -> maginput, maginputFieldNames, maginputUnits, info
% omniData -> minutely -> maginput, maginputFieldNames, maginputUnits, info
% matFilePath -> A string array recording the location where data is stored
%--------------------------------------------------------------------------
% Modified: 13th Feb 2017 
% Created : 6th Feb 2017
% Author  : Nithin Sivadas
% Ref     : 
% Notes   : 
%--------------------------------------------------------------------------

if nargin < 2
    dataStoreDir = [initialize_root_path,'LargeFiles\'];
end

if nargin < 1
    dateStr = '26 Mar 2008';
end
date = datenum(dateStr);
omniHourlyDir = [dataStoreDir,datestr(date,'yyyymmdd'),'\omni\hourly\'];
omni1minDir = [dataStoreDir,datestr(date,'yyyymmdd'),'\omni\1min\'];

% If the data doesn't exist in the appropriate folder, download them
if exist(omniHourlyDir)~=7||exist(omni1minDir)~=7
    download_omni_cdf(dateStr,dataStoreDir);
end

%% Generating Hourly Data
fileHourlyStr= get_files_in_folder(omniHourlyDir,'*.cdf');
VarNamesHourly1=[{'KP'};{'DST'}];
VarNamesHourly2=[{'N'};{'V'};{'Pressure'};{'BY_GSM'};{'BZ_GSM'}];

omniData.hourly.time = double...
    ([spdfcdfread([omniHourlyDir,fileHourlyStr{1}],'Variables',{'Epoch'},'ConvertEpochToDatenum', true)'...
    spdfcdfread([omniHourlyDir,fileHourlyStr{2}],'Variables',{'Epoch'},'ConvertEpochToDatenum', true)']);

omniData.hourly.maginput(:,1:2) =...
    double([cell2mat((spdfcdfread([omniHourlyDir,fileHourlyStr{1}],'Variables',VarNamesHourly1)));...
    cell2mat((spdfcdfread([omniHourlyDir,fileHourlyStr{2}],'Variables',VarNamesHourly1)))]);
omniData.hourly.maginput(:,3:7) =...
    double([cell2mat((spdfcdfread([omniHourlyDir,fileHourlyStr{1}],'Variables',VarNamesHourly2)));...
    cell2mat((spdfcdfread([omniHourlyDir,fileHourlyStr{2}],'Variables',VarNamesHourly2)))]);

info = spdfcdfinfo([omniHourlyDir,fileHourlyStr{1}]);
[tf, loc] = ismember(info.VariableAttributes.UNITS(:,1),[VarNamesHourly1',VarNamesHourly2']);
[~, p] = sort(loc(tf));
VarNamesHourlyInfoID = find(tf);
VarNamesHourlyInfoID = VarNamesHourlyInfoID(p);
omniData.hourly.maginputUnits(1,:) = [{'KP'};{'a.u.'}];
omniData.hourly.maginputUnits(2:7,:) = info.VariableAttributes.UNITS(VarNamesHourlyInfoID,:);

[tf, loc] = ismember(info.VariableAttributes.FIELDNAM(:,1),[VarNamesHourly1',VarNamesHourly2']);
[~, p] = sort(loc(tf));
VarNamesHourlyInfoID = find(tf);
VarNamesHourlyInfoID = VarNamesHourlyInfoID(p);
omniData.hourly.maginputFieldNames = info.VariableAttributes.FIELDNAM(VarNamesHourlyInfoID,:);

%% Calculating additional maginput variables for IRBEM [Incomplete]
 programDir = [initialize_root_path,...
        'energy-height-conversion\Tools\External Tools\Tsyganenko_Parameters\MagParameterProgram-rsw\'];
omniASCDir = [initialize_root_path,'LargeFiles\omni\ASC\'];
GW = get_tsyganenko_GW(str2num(datestr(date,'YYYY')),programDir,omniASCDir);  % This is minute-wise

%% Generating minute wise data
file1minStr= get_files_in_folder(omni1minDir,'*.cdf');
VarNames1min1=[{'SYM_H'}]; % Assume SYM_H = DST
VarNames1min2=[{'proton_density'};{'flow_speed'};{'Pressure'};{'BY_GSM'};{'BZ_GSM'}];

omniData.minutely.time = double...
    (spdfcdfread([omni1minDir,file1minStr{2}],'Variables',{'Epoch'},...
    'ConvertEpochToDatenum', true)');% Time
timeMin = omniData.minutely.time(1);
timeMax = omniData.minutely.time(end);
tempIndex = 1:1:length(GW.time);
timeIndex = crop_time(tempIndex,GW.time,timeMin,timeMax);
omniData.minutely.maginput(:,1) = interp1(omniData.hourly.time,...
    omniData.hourly.maginput(:,1),omniData.minutely.time,'nearest'); % KP
omniData.minutely.maginput(:,2) =...
    double(spdfcdfread([omni1minDir,file1minStr{2}],...
    'Variables',VarNames1min1)); %DST (SYM_H)
omniData.minutely.maginput(:,3:7) =...
    double(cell2mat((spdfcdfread([omni1minDir,file1minStr{2}],...
    'Variables',VarNames1min2)))); %N, V, P, BY_GSM, BZ_GSM
omniData.minutely.maginput(:,8:10) = GW.G(timeIndex,:); %G1,G2,G3
omniData.minutely.maginput(:,11:16) = GW.W(timeIndex,:); %W1-W6
omniData.minutely.maginput(:,17)  =...
    double(spdfcdfread([omni1minDir,file1minStr{2}],...
    'Variables',{'AL_INDEX'})); %AL index
omniData.minutely.maginput(:,18:25) = zeros(length(timeIndex),8); % For future use

info = spdfcdfinfo([omni1minDir,file1minStr{2}]);
[tf, loc] = ismember(info.VariableAttributes.UNITS(:,1),[VarNames1min1',VarNames1min2']);
[~, p] = sort(loc(tf));
VarNames1minInfoID = find(tf);
VarNames1minInfoID = VarNames1minInfoID(p);
omniData.minutely.maginputUnits(1,:) = [{'KP'};{'a.u.'}];
omniData.minutely.maginputUnits(2:7,:) = info.VariableAttributes.UNITS(VarNames1minInfoID,:);
omniData.minutely.maginputUnits(8:16,:) = [{'G1'},{'G2'},{'G3'},...
    {'W1'},{'W2'},{'W3'},{'W4'},{'W5'},{'W6'};...
    {'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'}]';
omniData.minutely.maginputUnits(17,:) = info.VariableAttributes.UNITS(29,:);
[tf, loc] = ismember(info.VariableAttributes.FIELDNAM(:,1),[VarNames1min1',VarNames1min2']);
[~, p] = sort(loc(tf));
VarNames1minInfoID = find(tf);
VarNames1minInfoID = VarNames1minInfoID(p);
omniData.minutely.maginputFieldNames(1,:) = [{'KP'};{'Kp*10 (3-h) index in minutes'}];
omniData.minutely.maginputFieldNames(2:7,:) = info.VariableAttributes.FIELDNAM(VarNames1minInfoID,:);
omniData.minutely.maginputUnits(8:16,:) = [{'G1'},{'G2'},{'G3'},...
    {'W1'},{'W2'},{'W3'},{'W4'},{'W5'},{'W6'};...
    {'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'},{'a.u.'}]';
omniData.minutely.maginputFieldNames(17,:) = info.VariableAttributes.FIELDNAM(41,:);
%%
omniData.minutely.info = ['Contains minute-wise data input for Tsyganenko models for ',datestr(date,'mmm YYYY')];
omniData.hourly.info = ['Contains hourly data input for Tsyganenko models for ',datestr(date,'YYYY')];
matFilePath = [omni1minDir,'omniData.mat'];
save(matFilePath,'omniData');

%% Creating a File List that records the the directories where the files are stored
processedFileList = [dataStoreDir,datestr(date,'yyyymmdd'),'\processedFileList.mat'];
if isfile(processedFileList)
    file=load(processedFileList);
    matlabFilePath = file.matlabFilePath;
    matlabFilePath.omni=matFilePath;
    save(processedFileList,'matlabFilePath','-append');
else
    matlabFilePath.omni=matFilePath;
    save(processedFileList,'matlabFilePath');
end

end
