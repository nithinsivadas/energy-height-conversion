function download_DASC_FITS(dateStr,storeDir)
%DOWNLOAD_DASC_FITS Downloads PokerFlat DASC FITS file

% Input
% dateStr   - 'dd mmm yyyy'; default : '26 Mar 2008'
% storeDir  - Directory where the FITS files should be downloaded to
%           - Default: '~\LargeFiles\DASC\'
% Ouput     - Downloads all FITS files not in the folder 

if nargin < 2
    storeDir = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep];
end

if nargin < 1
    dateStr = '26 Mar 2008';
end

host = 'optics.gi.alaska.edu';
remoteStoreLink = '/PKR/DASC/RAW/';
date = datenum(dateStr);
remoteFinalLink = [remoteStoreLink,datestr(date,'yyyy'),'/',datestr(date,'yyyymmdd')];
localStorePath = [storeDir,datestr(date,'yyyymmdd')];


dasc=ftp(host);
cd(dasc,remoteFinalLink);
remoteFileList = dir(dasc);
remoteFileListName = string(char(remoteFileList.name));

% Check if local Store Path exists
if isdir(localStorePath)
    localFileList = dir(localStorePath);
    localFileListName = string(char(localFileList.name));
    % Identifying identical files in remote and local directories
    index = arrayfun(@(k)sum(strcmp(remoteFileListName(k),localFileListName))>0, (1:length(remoteFileListName))','UniformOutput',false);
    indexMissingFiles = find(~cell2mat(index));
    % If all files are already downloaded then displaying a message
    if isempty(indexMissingFiles)
        disp('DASC data already downloaded');
    end
    % Downloading only the missing files
    for i=1:1:length(indexMissingFiles)
        mget(dasc,remoteFileListName(indexMissingFiles(i)),localStorePath);
    end
else
    mget(dasc,'*',localStorePath);
end


close(dasc);
end
