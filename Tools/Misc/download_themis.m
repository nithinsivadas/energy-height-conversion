function download_themis(dateStr,storeDir,probeStr,dataType)
%download_themis.m Downloads themis data
%--------------------------------------------------------------------------
% Input
%------
% dataStr - [dd mmm yyyy] A string indicating date for which you would like
%           to download the omni data
%           Default:'26 Mar 2008'
% storeDir- A string indicating the data storage directory in '...\github\'
%           It will store the data into ~\themis\state\
%           Default: '~\github\LargeFiles\'
% probeStr- A string listing the themis probes you would like the data from
%           Default: 'tha,thb,thc,thd,the' all of them
% dataType- A string listing the instrument from which data is to be downloaded 
%           'state', 'sst'  etc.
%           Default: 'state' - which is position of the satellite

%--------------------------------------------------------------------------
% Output
%-------
% Stores the cdf files in Store Dir
%--------------------------------------------------------------------------
% Modified: 12th Feb 2017 
% Created : 12th Feb 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------
    f = filesep;
    if nargin < 4
        dataType = 'state';
    end
    if nargin < 3
        probeStr = 'tha,thb,thc,thd,the';
    end
    if nargin < 2
        storeDir = [initialize_root_path,'LargeFiles\'];
    end
    if nargin < 1
        dateStr = '26 Mar 2008';
    end
    
    probe = strsplit(probeStr,',');
    date = datenum(dateStr);
    
    %% Downloading State of THEMIS
    if strcmp(dataType,'state')
        hWaitBar = waitbar(0,'Downloading themis data');
        host = 'spdf.gsfc.nasa.gov';
        thmFTP=ftp(host);
        for thisSC = 1:length(probe)
            thisProbe = char(probe(thisSC));
            custom_waitbar(hWaitBar,thisSC,length(probe),['Downloading ',...
                thisProbe,' : ',dataType,' ...'])
            dataStoreLink = ['/pub/data/themis/',(thisProbe),'/ssc/'];
            dataFinalLink = [dataStoreLink,datestr(date,'yyyy'),'/'];
            fileStr = [thisProbe,'_or_ssc_',datestr(date,'yyyymm'),'01_v01.cdf'];
            dataLocalDir = [storeDir,f,'themis',f,'state',f];
            if ~isfile([dataLocalDir,fileStr])
                cd(thmFTP,dataFinalLink);
                mget(thmFTP,fileStr,dataLocalDir);
            else 
                disp(['File ',fileStr,' is already present in ',dataLocalDir]);
            end
            close(thmFTP)
        end
        delete(hWaitBar);
    end

end

