%% Find processed AMISR Exp Data 
clear variables
% Step 1: Get all experiment folder names
expDatabaseFilePath = 'G:\Team Drives\Semeter-Research in Progress\All AMISR Experiments\RawData_Folders_Exps_Times_by_Radar_9_April_2018.h5';
webDatabase = 'G:\Team Drives\Semeter-Research in Progress\All AMISR Experiments\amisrWebDatabase_May_2018.h5';
superMagDatabase = 'G:\My Drive\Research\Projects\Paper 2\Data\Events List\20180605-20-14-substorms.txt';
%% Processing radar database
[data,h5DataLoc] = h5read_data(expDatabaseFilePath);
% expNames=regexprep(num2cell(data.ExpNames.Names.ExpName,2),'\s','');
% expMountPaths=regexprep(num2cell(data.MountPaths.Path.MountPath,2),'\s','');
expNames=deblank(num2cell(data.ExpNames.Names.ExpName,2));
expMountPaths=deblank(num2cell(data.MountPaths.Path.MountPath,2));

radarDatabase = data.Radars.PFISR;
radarDatabase=get_exp_folder(radarDatabase);
radarDatabase.nExpId=radarDatabase.nExpId+1;%correcting for matlab indices
radarDatabase.nMountPathId=radarDatabase.nMountPathId+1; %correcting for matlab indices

%% Input
interestedExpNames=expNames(find(strncmpi(expNames,'Sporadic',8) | ...
    strncmpi(expNames,'MSWind',6) | strncmpi(expNames,'Themis',6) | ...
    strncmpi(expNames,'Lyons',5)));
% interestedExpNames=expNames;

%% Results
selectedExperiments = get_experiment_path(interestedExpNames,expNames,expMountPaths,radarDatabase);
% [folderArray,mount]=create_file_list(selectedExperiments,expMountPaths,false,'\Setup\*.exp');

%% Check if processed or unprocessed
webData=h5read_data(webDatabase);
nInterestedExp = length(interestedExpNames);
for iExp = 1:1:nInterestedExp
    nExpID = length(selectedExperiments(iExp).folderName);
    for iExpID = 1:1:nExpID
        indx=find(strcmp(deblank(string(webData.PFISR.expId)),selectedExperiments(iExp).folderName(iExpID)));
        try
            selectedExperiments(iExp).status(iExpID,1)={(deblank(string(webData.PFISR.processingStatus(indx(1)))))};
            selectedExperiments(iExp).startTime(iExpID,1)={unix_to_matlab_time(webData.PFISR.startTime(indx(1)))};
            selectedExperiments(iExp).endTime(iExpID,1)={unix_to_matlab_time(webData.PFISR.endTime(indx(1)))};
        catch ME
            selectedExperiments(iExp).status(iExpID,1)={'Unknown'};
            selectedExperiments(iExp).startTime(iExpID,1)={0};
            selectedExperiments(iExp).endTime(iExpID,1)={0};
        end
    end
end

%% Print out
k = 1;
for iExp = 1:1:nInterestedExp
    nExpID = length(selectedExperiments(iExp).folderName);
    for iExpID = 1:1:nExpID
        writeTable(k,1) = {string(interestedExpNames(iExp))};
        writeTable(k,2) = {string(selectedExperiments(iExp).folderName(iExpID,1))}; 
        writeTable(k,3) = selectedExperiments(iExp).status(iExpID,1); 
        writeTable(k,4) = selectedExperiments(iExp).startTime(iExpID,1);
        writeTable(k,5) = selectedExperiments(iExp).endTime(iExpID,1);
        k=k+1;
    end
end
writeTable(find(cellfun(@isempty,writeTable(:,3))),3)={'Unknown'};
% processedIndx=find(strcmp(string(writeTable(:,3)),'Processed and Available - Calibrated'));
% T = cell2table(writeTable(:,1:3),'VariableNames',{'ExpName','ExpID','Status'});
% writetable(T,'Lyons.xlsx');

%% SuperMag Database
substorms = dlmread(superMagDatabase,'',68,0);
tempSubstormDate = datenum(datetime([substorms(:,1:5),ones(size(substorms,1),1)]));
substormList(:,1) = cellstr(datestr(datetime([substorms(:,1:5),ones(size(substorms,1),1)])));
substormList(:,2) = num2cell(substorms(:,6));
substormList(:,3) = num2cell(substorms(:,7));
%%
for i=1:1:length(substormList(:,1))
    try 
        substormList(i,4)={string(writeTable{cell2mat(writeTable(:,4))<=tempSubstormDate(i) & cell2mat(writeTable(:,5))>=tempSubstormDate(i),2})};
        substormList(i,5)={string(writeTable{cell2mat(writeTable(:,4))<=tempSubstormDate(i) & cell2mat(writeTable(:,5))>=tempSubstormDate(i),1})};
        substormList(i,6)={string(writeTable{cell2mat(writeTable(:,4))<=tempSubstormDate(i) & cell2mat(writeTable(:,5))>=tempSubstormDate(i),3})};
    catch ME
        substormList(i,4) = {"nan"};
    end
end
substormTable = cell2table(substormList,'VariableNames',{'Time','MLAT','MLT','ExpID','ExpName','processingStatus'});
writetable(substormTable,'substorm_2001_to_2017.txt');

%% Functoins
function [folderArray,mount]=create_file_list(selectedExperiments,expMountPaths,setSaveFileOn,sizeCalculationFilesWildcard)

if nargin<4
    sizeCalculationFilesWildcard='\*';
end

if nargin<3
    setSaveFileOn=true;
end

exp=struct2cell(selectedExperiments');
folderArray.mountPathID =[];
folderArray.name ={};
for i=1:1:length(exp)
    folderArray.mountPathID=[folderArray.mountPathID; exp{1,i}'];
    folderArray.name=[folderArray.name; exp{2,i}];
end
uniqueMountPaths = unique(folderArray.mountPathID);
for iMount = 1:1:length(uniqueMountPaths)
    indx = find(folderArray.mountPathID==uniqueMountPaths(iMount));
    mount(iMount).folder = folderArray.name(indx);
    mount(iMount).path = expMountPaths{uniqueMountPaths(iMount)};
    mount(iMount).fileListName = extractBetween(mount(iMount).path,'/Volumes/','/DataAMISRPoker');
    if setSaveFileOn
        fileID = fopen(string(strcat(mount(iMount).fileListName,'.txt')),'w');
        fprintf(fileID,'%s\n',string(mount(iMount).folder));
        fclose(fileID);
    end
    mount(iMount).GBytes = calculate_size(mount(iMount).folder,uniqueMountPaths(iMount),sizeCalculationFilesWildcard); 
end

end

function fileName=get_experiment_path(expNameCell,expNames,mountPaths,radarDatabase)

for iExp=1:1:length(expNameCell)
    expId(iExp)=find(strcmp(string(expNames),expNameCell(iExp)));
    indx=radarDatabase.nExpId==expId(iExp);
    thisExp(iExp).nMountPathId=radarDatabase.nMountPathId(indx);
    thisExp(iExp).folderName = radarDatabase.expFolderName(indx);
end
fileName=thisExp;
end
function radarDatabase=get_exp_folder(radarDatabase)

radarDatabase.expFolderName =...
    compose('%04i%02i%02i.%03i',...
    [radarDatabase.nyear',radarDatabase.nmonth',radarDatabase.nday',radarDatabase.nset']);

end

function [data,h5DataLoc]=h5read_data(h5FilePath)
info=h5info(h5FilePath);
nGroups=length(info.Groups);
for thisGroup=1:1:nGroups
    for thisDataSet = 1:1:length(info.Groups(thisGroup).Datasets)
    
    thisGroupName=info.Groups(thisGroup).Name(2:end);
    thisDataSetName=info.Groups(thisGroup).Datasets(thisDataSet).Name;
    
    h5DataLoc{thisGroup,thisDataSet} = ...
        strcat('/',thisGroupName,'/',thisDataSetName);
    
    temp=h5read(h5FilePath,h5DataLoc{thisGroup,thisDataSet});
    if isstruct(temp)
        data.(thisGroupName).(thisDataSetName) = ...
            structfun(@transpose,...
            temp,...
            'UniformOutput',false);
    else
        data.(thisGroupName).(thisDataSetName) = ...
            transpose(temp);
    end
    end
end
end

function GBytes=calculate_size(folderNameArray,mountPathID,wildCard)
if(nargin<3)
    wildCard = '\*';
end
GBytes=zeros(length(folderNameArray),1);
%mountPathID - just one integer
serverPath = NAS_server(mountPathID);

    for i=1:1:length(folderNameArray)
        thisExpFolderStr=folderNameArray{i};
        myfolderinfo = dir([serverPath,thisExpFolderStr,wildCard]);
        myfolderinfoCell=(struct2cell(myfolderinfo)');
        GBytes(i)=sum(cell2mat(myfolderinfoCell(:,4)))*10^-9; %GB
    end

end
function serverPath=NAS_server(mountPathID)
    mountPathID
    switch mountPathID
        case 2
        serverPath = '\\128.18.144.181\pfisr_001\Data AMISR Poker\';
        case 3
        serverPath = '\\128.18.144.182\pfisr_002\Data AMISR Poker\'; 
        case 4
        serverPath = '\\128.18.144.187\pfisr_005\Data AMISR Poker\';
        case 5
        serverPath = '\\128.18.144.192\pfisr_006\Data AMISR Poker\';         
    end
end