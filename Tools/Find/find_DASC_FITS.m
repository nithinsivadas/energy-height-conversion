function [timeStamp, wavelength, error]=find_DASC_FITS(timeStr)
%FIND_DASC_FITS Checks if there is data during the day of timeStr in UT. 

% Input
% timeStr    - 'dd mmm yyyy'; default : '26 Mar 2008'

% Ouput     
% timeStamp  - All time instances where data is available in [matlab time]
% wavelength - Corresponding wavelength of data, in [nano meters] 

% Created on: 13th Aug 2019
% Created by: Nithin Sivadas, Boston University

if nargin < 1
    timeStr = '26 Mar 2008';
end

host = 'optics.gi.alaska.edu';
remoteStoreLink = '/PKR/DASC/RAW/';
date = datenum(timeStr);
remoteFinalLink = [remoteStoreLink,datestr(date,'yyyy'),'/',datestr(date,'yyyymmdd')];

dasc=ftp(host);
ME = [];
try
    cd(dasc,remoteFinalLink);
    remoteFileList = dir(dasc);
    remoteFileListName = string(char(remoteFileList.name));

    try
        splitFileNameStr = cellfun(@strsplit,cellstr(remoteFileListName),repmat(cellstr(filesep),length(remoteFileListName),1));
        [timeStamp,wavelength] = fitsfiletimestamp(splitFileNameStr);

    catch ME
        timeStamp = nan;
        wavelength = nan;
    end
catch ME
    timeStamp = nan;
    wavelength = nan;
end   
error = ME;
close(dasc);
end