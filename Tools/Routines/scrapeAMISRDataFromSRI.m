%% Scrape AMISR Experiment Data
clear variables
% Step 1: Get all experiment folder names
expDatabaseFilePath = 'G:\Team Drives\Semeter-Research in Progress\All AMISR Experiments\RawData_Folders_Exps_Times_by_Radar_9_April_2018.h5';

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
interestedExpNames=expNames(find(strncmpi(expNames,'MSWinds23',9)));
% interestedExpNames=expNames;

%% Results
selectedExperiments = get_experiment_path(interestedExpNames,expNames,expMountPaths,radarDatabase);
[folderArray,mount]=create_file_list(selectedExperiments,expMountPaths,false,'\Setup\*.exp');

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
    
    data.(thisGroupName).(thisDataSetName) = ...
        structfun(@transpose,...
        h5read(h5FilePath,h5DataLoc{thisGroup,thisDataSet}),...
        'UniformOutput',false);
    
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