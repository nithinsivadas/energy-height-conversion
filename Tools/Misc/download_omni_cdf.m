function download_omni_cdf(dateStr,storeDir,dataType)
%download_omni_cdf.m Downloads hourly and minute-wise omni values
%--------------------------------------------------------------------------
% Input
%------
% dataStr - [dd mmm yyyy] A string indicating date for which you would like
%           to download the omni data
%           Default:'26 Mar 2008'
% storeDir- A string indicating the data storage directory in '...\github\'
%           It will store the hourly data in ~\yyyymmdd\omni\hourly\
%           and 1 min data in ~yyyymmdd\omni\1min\
%           Default: '~\github\LargeFiles\'
%--------------------------------------------------------------------------
% Output
%-------
% Stores the cdf files in Store Dir
%--------------------------------------------------------------------------
% Modified: 6th Feb 2017 
% Created : 6th Feb 2017
% Author  : Nithin Sivadas
% Ref     : 
% Notes   : The CDF file is stored in an awkard location - under each event
%           date. Only processed .mat files need to be stored here. Change
%           this.
%--------------------------------------------------------------------------
if nargin < 3
    dataType = 'cdf';
end
if nargin < 2
    storeDir = [initialize_root_path,'LargeFiles\'];
end
if nargin < 1
    dateStr = '26 Mar 2008';
end
date = datenum(dateStr);

if strcmp(dataType,'cdf')
    host = 'spdf.gsfc.nasa.gov';
    dataStoreLink = '/pub/data/omni/omni_cdaweb/';

    dataFinalLinkHourly = [dataStoreLink,'hourly/',datestr(date,'yyyy'),'/'];
    dataFinalLink1min = [dataStoreLink,'hro_1min/',datestr(date,'yyyy'),'/'];
    omni2=ftp(host);
    cd(omni2,dataFinalLinkHourly);
    mget(omni2,'*.cdf',[storeDir,datestr(date,'yyyymmdd'),'\omni\hourly\']);
    cd(omni2,dataFinalLink1min);
    if str2num(datestr(date,'mm'))~=1
    prevmm=datestr([datestr(date,'yyyy'),'/',...
        num2str(str2num(datestr(date,'mm'))-1),'/',num2str(01)],'mm');
    mget(omni2,['omni_hro_1min_',datestr(date,'yyyy'),prevmm,'01_v01.cdf'],...
        [storeDir,datestr(date,'yyyymmdd'),'\omni\1min\']); % previous month
    else
    prevyy=num2str((str2num(datestr(date,'yyyy'))-1));
    dataPrevLink1min = [dataStoreLink,'hro_1min/',prevyy,'/'];
    cd(omni2,dataPrevLink1min);
    mget(omni2,['omni_hro_1min_',prevyy,num2str(12),'01_v01.cdf'],...
        [storeDir,datestr(date,'yyyymmdd'),'\omni\1min\']); % previous month
    cd(omni2,dataFinalLink1min);
    end     
    mget(omni2,['omni_hro_1min_',datestr(date,'yyyymm'),'01_v01.cdf'],...
        [storeDir,datestr(date,'yyyymmdd'),'\omni\1min\']); % current month
    close(omni2);
end

if strcmp(dataType,'asc')
    host = 'nssdcftp.gsfc.nasa.gov';
    dataStoreLink = '/pub/data/omni/high_res_omni/';

    omni2=ftp(host);
    cd(omni2,dataStoreLink);
    prevyy=num2str((str2num(datestr(date,'yyyy'))-1));
    data1minPrevStr=['omni_min',prevyy,'.asc'];
    data1minStr=['omni_min',datestr(date,'yyyy'),'.asc'];
    mget(omni2,data1minPrevStr,[storeDir,'\omni\ASC\']);
    mget(omni2,data1minStr,[storeDir,'\omni\ASC\']);
    
    data5minPrevStr=['omni_5min',prevyy,'.asc'];
    data5minStr=['omni_5min',datestr(date,'yyyy'),'.asc'];
    mget(omni2,data5minPrevStr,[storeDir,'\omni\ASC\']);
    mget(omni2,data5minStr,[storeDir,'\omni\ASC\']);
    
    % Downloading CDF files that will let you create kp and dst files
    dataFinalLinkHourly = [dataStoreLink,'../omni_cdaweb/hourly/',datestr(date,'yyyy'),'/'];
    cd(omni2,dataFinalLinkHourly);
    mget(omni2,'*.cdf',[storeDir,'\omni\CDF\']);
    dataFinalLinkHourly = [dataStoreLink,'../omni_cdaweb/hourly/',prevyy,'/'];
    cd(omni2,dataFinalLinkHourly);
    mget(omni2,'*.cdf',[storeDir,'\omni\CDF\']);
    
    close(omni2);
end


end

