%% Routine to calculate the size of experiments on SRI Servers
% clear variables;
clear TBytes GBytes
% dateMinStr = '1 Jan 2007';
% dateMaxStr = '31 March 2016';
dateMinStr = '1 May 2016';
dateMaxStr = '1 April 2018';
% wildCard = '\*.dt0.h5';
wildCard = '\*';

[table,tableHeader] = get_amisr_experiment_times(dateMinStr,dateMaxStr);
DregionExps = {'MSWinds01';'MSWinds02';'MSWinds03';'MSWinds04';'MSWinds05';...
    'MSWinds15';'MSWinds16';'MSWinds17';'MSWinds17_3dt';'MSWinds18';'MSWinds19';...
    'MSWinds21';'MSWinds22';'MSWinds23';'MSWinds23_3dt';'MSWinds23hr';...
    'MSWinds23m';'MSWinds23_dt013';'MSWinds26.v01';'MSWinds26.v02';'MSWinds26.v03';...
    'MSWinds24';'MSWinds25';'MSWindsrand_v01';'MSWindsrand_v02'};
conjunctionExps = {'Themis30';'Themis31';'Themis36';'Lyons02';'Lyons13';'Lyons30_3dt'};
expName={'MSWinds23'};
TBytes=zeros(length(expName),1);
for j=1:1:length(expName)
    expIndx = find((strcmp(table(:,2),expName(j))));
    expFolderStr = table(expIndx,1);
    expFolderList(j)={table(expIndx,:)};
    GBytes = zeros(length(expIndx),1);
    for i=1:1:length(expIndx)
        thisExpFolderStr=expFolderStr{i};
        serverPath = NAS_server(thisExpFolderStr);
        myfolderinfo = dir([serverPath,thisExpFolderStr,wildCard]);
        myfolderinfoCell=(struct2cell(myfolderinfo)');
        GBytes(i)=sum(cell2mat(myfolderinfoCell(:,4)))*10^-9; %GB
    end
    TBytes(j) = sum(GBytes)*10^-3;
    disp(['Exp: ',cellstr(expName(j)),' ',dateMinStr,' to ',dateMaxStr,' Size:',num2str((TBytes(j)),3),' TB ','Files:',wildCard]);
end

function serverPath=NAS_server(thisExpFolderStr)
    thisExpDate=datenum(thisExpFolderStr(1:8),'yyyymmdd');
    if(thisExpDate <= datenum('20141130','yyyymmdd'))
        serverPath = '\\128.18.144.181\pfisr_001\Data AMISR Poker\';
    elseif(thisExpDate <= datenum('20160111','yyyymmdd'))
        serverPath = '\\128.18.144.182\pfisr_002\Data AMISR Poker\'; 
    elseif(thisExpDate <= datenum('20171031','yyyymmdd'))
        serverPath = '\\128.18.144.187\pfisr_005\Data AMISR Poker\';
    else
        serverPath = '\\128.18.144.192\pfisr_006\Data AMISR Poker\';    
    end
end

