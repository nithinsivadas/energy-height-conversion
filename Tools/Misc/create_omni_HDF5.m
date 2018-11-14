function [data,status] = create_omni_HDF5(h5FileStr,localStorePath,setCalculateGW)
%create_omni_HDF5
%
if nargin<3
    setCalculateGW = true;
end

if nargin < 2
    localStorePath = [initialize_root_path,'LargeFiles',filesep,'omni',filesep,'ASC'];
end

if nargin <1
    h5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';
end

% Setting up remote path 
remoteLinkOmniHighRes = '/pub/data/omni/high_res_omni';
remoteLinkOmni2 = '/pub/data/omni/low_res_omni/extended';

%% OMNI High_res_download
download_omni_high_res(remoteLinkOmniHighRes, localStorePath);

%% OMNI2
download_omni_2(remoteLinkOmni2, localStorePath);
status{1,1} = 'Download complete';

%% Load ASCII Filess
% OMNI High_Res
create_omni_HDF5_file(localStorePath,h5FileStr,setCalculateGW);
status{2,1} = 'Loaded ASCII Files';

end
function create_omni_HDF5_file(localStorePath, h5FileStr, setCalculateGW)
    
    if nargin<3
        setCalculateGW = true;
    end
    
    files=dir(localStorePath);
    fileListCell = struct2cell(files);
    formatSpecOmni1 = '%4d %4d %3d %3d %3d %3d %4d %4d %4d %7d %7d %6.2f %7d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.1f %8.1f %8.1f %8.1f %7.2f %9.0f %6.2f %7.2f %7.2f %6.1f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %6d %6d %6d %6d %6d %6d %6d %7.2f %5.1f';
    formatSpecOmni2 = '%4d %4d %3d %5d %3d %3d %4d %4d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %9.0f %6.1f %6.0f %6.1f %6.1f %6.3f %6.2f %9.0f %6.1f %6.0f %6.1f %6.1f %6.3f %7.2f %7.2f %6.1f %3d %4d %6d %5d %10.2f %9.2f %9.2f %9.2f %9.2f %9.2f %3d %4d %6.1f %6.1f %6d %6d %5.1f %6.2f';
    fileListMinStr=strcat(fileListCell(2,strncmp(fileListCell(1,:)','omni_min',8))',filesep,fileListCell(1,strncmp(fileListCell(1,:)','omni_min',8))');
    n = length(fileListMinStr);
    k=1;
    l=1;
    multiWaitbar('Stitch omni_high_res data',0);
    
    h5Headers = {'/Time','/Timeshift',...
        '/BField/BxGSE', '/BField/ByGSE','/BField/BzGSE',...
        '/BField/ByGSM', '/BField/BzGSM',...
        '/Velocity/V','/Velocity/VxGSE','/Velocity/VyGSE','/Velocity/VzGSE',...
        '/ProtonDensity','/Temperature','/FlowPressure','/EField','/PlasmaBeta','/AlfvenMach','/AlphaProtonRatio',...
        '/Indices/AE','/Indices/AL','/Indices/AU',...
        '/Indices/SYM_D','/Indices/SYM_H','/Indices/Kp',...
        '/Error/sigmaTemperature','/Error/sigmaProtonDensity',...
        '/Error/sigmaVelocity','/Error/sigmaAlphaProtonRatio',...
        '/TSY/G1','/TSY/G2','/TSY/G3',...
        '/TSY/W1','/TSY/W2','/TSY/W3','/TSY/W4','/TSY/W5','/TSY/W6'};
    
    h5Description = {'UNIX Date-Time','Timeshift from L1 to Stagnation Point',...
        'Bx - same for GSE and GSM', 'By - GSE','Bz - GSE',...
        'By - GSM', 'Bz - GSM',...
        'Flow speed','Velocity along xGSE','Velocity along yGSE','Velocity along zGSE',...
        'Proton density','Temperature','Flow pressure','Electric field','Plasma beta','Alfven Mach Number','Alpha/Proton Ratio (1hr)',...
        'AE Index','AL Index','AU Index',...
        'SYM_D Index','SYM_H Index','Planetary Geomagnetic Activity Index (e.g. 3+ = 33, 6- = 57, 4 = 40) [1hr]',...
        'sigmaT [1hr]','sigmaN [1hr]',...
        'sigmaV [1hr]','sigma_Na_Np [1hr]',...
        'TSY Parameter G1','TSY Parameter G2','TSY Parameter G3',...
        'TSY Parameter W1','TSY Parameter W2','TSY Parameter W3','TSY Parameter W4','TSY Parameter W5','TSY Parameter W6'};
    
    h5Units = {'UNIX Time [days]','[sec]',...
        '[nT]', '[nT]','[nT]',...
        '[nT]', '[nT]',...
        '[km/s]','[km/s]','[km/s]','[km/s]',...
        '[n/cc]','[k]','[nPa]','[mV/m]','[a.u.]','[a.u.]','[a.u.]',...
        '[nT]','[nT]','[nT]',...
        '[nT]','[nT]','[a.u.]',...
        '[deg k]','[n/cc]',...
        '[km/s]','[a.u.]',...
        '[a.u.]','[a.u.]','[a.u.]',...
        '[a.u.]','[a.u.]','[a.u.]','[a.u.]','[a.u.]','[a.u.]'};
    
    q = length(h5Headers);
    
   
    for j = 1:1:q
            h5create(h5FileStr,h5Headers{j},[1 Inf],'ChunkSize',[1 100],'Deflate',9);
            h5writeatt(h5FileStr,h5Headers{j},'Descriptions',h5Description{j});
            h5writeatt(h5FileStr,h5Headers{j},'Units',h5Units{j});
    end
    
    for i=1:1:n
        multiWaitbar('Stitch omni_high_res data','Value',i/n);
        yr = str2num(fileListMinStr{i}(regexp(fileListMinStr{i},'_min')+4:regexp(fileListMinStr{i},'.asc')-1));
        
        if yr >= 1995
        
        %% omni1 
        omni1 = load_ascii_files(fileListMinStr{i},formatSpecOmni1);
        
        kwidth = length(omni1{1,1});
        kend = k+kwidth-1;
        
        time1 = datetime(omni1{1,1},ones(kwidth,1),...
            omni1{1,2},omni1{1,3},omni1{1,4},...
            zeros(kwidth,1));
        
        %% omni2
        fileOmni2Str=strcat(fileListCell(2,strncmp(fileListCell(1,:)',['omni2_',num2str(yr)],10))',filesep,fileListCell(1,strncmp(fileListCell(1,:)',['omni2_',num2str(yr)],10))');
        omni2 = load_ascii_files(fileOmni2Str{1},formatSpecOmni2);     
        lwidth = length(omni2{1,1});
        time2 = datetime(omni2{1,1},ones(lwidth,1),...
            omni2{1,2},omni2{1,3},zeros(lwidth,1),...
            zeros(lwidth,1));
        
        data = select_required_parameters(omni1, omni2, time1, time2, yr, setCalculateGW);
                
        for j = 1:1:q
                h5write(h5FileStr,h5Headers{j},data(:,j)',[1 k],[1 kwidth]);
        end 
        
        k = kend+1;
        end
        
    end
    multiWaitbar('Stitch omni_high_res data',1);
    
end

function data=select_required_parameters(omni1, omni2, time1, time2, yyyy, setCalculateGW)
        
        if nargin<6
            setCalculateGW = true;
        end
        
        ptime2 = posixtime(time2);
        data(:,1) = posixtime(time1);
        data(:,2) = omni1{1,10}; data(data(:,2)==999999,2)=nan; %Timeshift
        data(:,3) = omni1{1,15}; data(data(:,3)==9999.99,3)=nan; %Bx GSE
        data(:,4) = omni1{1,16}; data(data(:,4)==9999.99,4)=nan; %By GSE
        data(:,5) = omni1{1,17}; data(data(:,5)==9999.99,5)=nan; %Bz GSE
        data(:,6) = omni1{1,18}; data(data(:,6)==9999.99,6)=nan; %By GSM
        data(:,7) = omni1{1,19}; data(data(:,7)==9999.99,7)=nan; %By GSM
        data(:,8) = omni1{1,22}; data(data(:,8)==99999.9,8)=nan; %Flow Speed
        data(:,9) = omni1{1,23}; data(data(:,9)==99999.9,9)=nan; %Vx
        data(:,10) = omni1{1,24};data(data(:,10)==99999.9,10)=nan; %Vy
        data(:,11) = omni1{1,25};data(data(:,11)==99999.9,11)=nan; %Vz
        data(:,12) = omni1{1,26};data(data(:,12)==999.99,12)=nan; %Proton Density
        data(:,13) = omni1{1,27};data(data(:,13)==9999999.0,13)=nan; %Temperature
        data(:,14) = omni1{1,28};data(data(:,14)==99.99,14)=nan; %Flow pressure
        data(:,15) = omni1{1,29};data(data(:,15)==999.99,15)=nan; %Electric field
        data(:,16) = omni1{1,30};data(data(:,16)==999.99,16)=nan; % Plasma Beta
        data(:,17) = omni1{1,31};data(data(:,17)==999.9,17)=nan; % Alfven Mach Number
        data(:,18) = interp1(ptime2,omni2{1,28},data(:,1),'nearest'); data(data(:,18)==9.999,18)=nan; % Alpha Proton Ratio
        data(:,19) = omni1{1,38};data(data(:,19)==9999.99,19)=nan; % AE-Index
        data(:,20) = omni1{1,39};data(data(:,20)==9999.99,20)=nan; % AL-Index
        data(:,21) = omni1{1,40};data(data(:,21)==9999.99,21)=nan; % AU-Index
        data(:,22) = omni1{1,41};data(data(:,22)==9999.99,22)=nan; % SYM/D-Index
        data(:,23) = omni1{1,42};data(data(:,23)==9999.99,23)=nan; % SYM/H-Index
        data(:,24) = interp1(ptime2,double(omni2{1,39}),data(:,1),'nearest');data(data(:,24)==99.0,24)=nan; % Kp Index
        data(:,25) = interp1(ptime2,omni2{1,30},data(:,1),'nearest'); data(data(:,25)==9999999.0,25)=nan; % Sigma T
        data(:,26) = interp1(ptime2,omni2{1,31},data(:,1),'nearest'); data(data(:,26)==999.9,26)=nan; % Sigma N
        data(:,27) = interp1(ptime2,omni2{1,32},data(:,1),'nearest'); data(data(:,27)==9999.0,27)=nan; % Sigma V
        data(:,28) = interp1(ptime2,omni2{1,33},data(:,1),'nearest'); data(data(:,28)==999.9,28)=nan; % Sigma Alpha Proton Ratio
        
        if setCalculateGW==true
            GW = get_tsyganenko_GW_1(yyyy);
            data(:,29:31) = interp1(GW.time, GW.G, unixtime2matlab(data(:,1)),'nearest');
            data(:,32:37) = interp1(GW.time, GW.W, unixtime2matlab(data(:,1)),'nearest');
        else
            GW.time = data(:,1);
            data(:,29:31) = nan(length(GW.time),3);
            data(:,32:37) = nan(length(GW.time),6);
        end
        

end

function data=load_ascii_files(loadFile, format)
fileID = fopen(loadFile,'r');
data = textscan(fileID, format);
fclose(fileID);
end

function download_omni_high_res(remoteLink, localStorePath)

host = 'spdf.gsfc.nasa.gov';
remoteFinalLink = remoteLink;
omni=ftp(host);
cd(omni,remoteFinalLink);
remoteFileList = dir(omni);
remoteFileListCell = struct2cell(remoteFileList);
remoteFileListNameIndx = strncmp(remoteFileListCell(1,:),'omni_min',8);
remoteFileListName = string(remoteFileListCell(1,remoteFileListNameIndx));

% Download 1-min files
if isdir(localStorePath)
    localFileList = dir(localStorePath);
    localFileListName = strtrim(string(char(localFileList.name)));
    % Identifying identical files in remote and local directories
    index = arrayfun(@(k)sum(strcmp(remoteFileListName(k),localFileListName))>0, (1:length(remoteFileListName))','UniformOutput',false);
    indexMissingFiles = find(~cell2mat(index));
    % Making sure current year data is always downloaded
    currentYearIndx = find((strncmp(remoteFileListName,strcat('omni_min',string(year(datetime))),12)));
    if sum(indexMissingFiles==currentYearIndx)==0
        indexMissingFiles(end+1)=currentYearIndx;
    end
    % If all files are already downloaded then displaying a message
    if isempty(indexMissingFiles)
        disp('DASC data already downloaded');
    end
    %Downloading the missing files
    urlFile = 'omni_1min_URLs.txt';
    urlFilePath = [localStorePath,filesep,urlFile];
    urls = strcat('ftp://',host,remoteFinalLink,'/',remoteFileListName(indexMissingFiles));
    
    fileID = fopen(urlFilePath,'w'); fprintf(fileID,'%s\r\n',urls');fclose(fileID);
    if isunix
    [status,cmdout]=unix(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,' -i ',urlFilePath]);
    else
    [status,cmdout]=system(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,' -i ',urlFilePath]);
    end
else
    mkdir(localStorePath);
    urlFile = 'omni_1min_URLs.txt';
    urlFilePath = [localStorePath,filesep,urlFile];
    urls = strcat('ftp://',host,remoteFinalLink,'/',remoteFileListName);
    fileID = fopen(urlFilePath,'w'); fprintf(fileID,'%s\r\n',urls');fclose(fileID);
    if isunix
        [status,cmdout]=...
            unix(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,'-i ',urlFilePath]);
    else
        [status,cmdout]=...
            system(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,'-i ',urlFilePath]);
    end    
end

close(omni);

end

function download_omni_2(remoteLink, localStorePath)

%% OMNI-2 Download 
host = 'spdf.gsfc.nasa.gov';
omni=ftp(host);
remoteFinalLinkHourly = remoteLink;
cd(omni,remoteFinalLinkHourly);
remoteFileList = dir(omni);
remoteFileListCell = struct2cell(remoteFileList);
remoteFileListNameIndx = strncmp(remoteFileListCell(1,:),'omni2_',6);
remoteFileListNameIndx = remoteFileListNameIndx&~endsWith(remoteFileListCell(1,:),'all_years.dat');
remoteFileListName = string(remoteFileListCell(1,remoteFileListNameIndx));

% Download 1-hour files
if isdir(localStorePath)
    localFileList = dir(localStorePath);
    localFileListName = strtrim(string(char(localFileList.name)));
    % Identifying identical files in remote and local directories
    index = arrayfun(@(k)sum(strcmp(remoteFileListName(k),localFileListName))>0, (1:length(remoteFileListName))','UniformOutput',false);
    indexMissingFiles = find(~cell2mat(index));
    % Making sure current year data is always downloaded
    currentYearIndx = find((strncmp(remoteFileListName,strcat('omni2_',string(year(datetime))),10)));
    % Remove current & previous year?
    if sum(indexMissingFiles==currentYearIndx)==0
        indexMissingFiles(end+1)=currentYearIndx;
    end
    % If all files are already downloaded then displaying a message
    if isempty(indexMissingFiles)
        disp('DASC data already downloaded');
    end
    %Downloading the missing files
    urlFile = 'omni2_URLs.txt';
    urlFilePath = [localStorePath,filesep,urlFile];
    urls = strcat('ftp://',host,remoteFinalLinkHourly,'/',remoteFileListName(indexMissingFiles));
    fileID = fopen(urlFilePath,'w'); fprintf(fileID,'%s\r\n',urls');fclose(fileID);
    if isunix
    [status,cmdout]=unix(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,' -i ',urlFilePath]);
    else
    [status,cmdout]=system(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,' -i ',urlFilePath]);
    end
else
    mkdir(localStorePath);
    urlFile = 'omni2_URLs.txt';
    urlFilePath = [localStorePath,filesep,urlFile];
    urls = strcat('ftp://',host,remoteFinalLinkHourly,'/',remoteFileListName);
    fileID = fopen(urlFilePath,'w'); fprintf(fileID,'%s\r\n',urls');fclose(fileID);
    if isunix
        [status,cmdout]=...
            unix(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,'-i ',urlFilePath]);
    else
        [status,cmdout]=...
            system(['aria2c --allow-overwrite true -V -c -j 50 ','-d ',localStorePath,'-i ',urlFilePath]);
    end    
end
close(omni);
end
