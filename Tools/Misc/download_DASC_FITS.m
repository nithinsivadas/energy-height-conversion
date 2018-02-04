function download_DASC_FITS(dateStr,storeDir)
%DOWNLOAD_DASC_FITS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    storeDir = [initialize_root_path,'LargeFiles\DASC\'];
end

if nargin < 1
    dateStr = '26 Mar 2008';
end

host = 'optics.gi.alaska.edu';
dataStoreLink = '/PKR/DASC/RAW/';
date = datenum(dateStr);
dataFinalLink = [dataStoreLink,datestr(date,'yyyy'),'/',datestr(date,'yyyymmdd')];

dasc=ftp(host);
cd(dasc,dataFinalLink);
mget(dasc,'*',[storeDir,datestr(date,'yyyymmdd')]);
close(dasc);
end



