function [data, error] = get_DASC_FITS_times(timeMinStr, timeMaxStr)
%GET_DASC_FITS_TIMES Stores details of images available on DASC server. 

% Input
% timeMinStr    - 'dd mmm yyyy'; default : '31 Dec 2017'
% timeMaxStr    - 'dd mmm yyyy'; default : '19 Jan 2018'

% Ouput     
% data
% --> time       - All time instances where data is available in [matlab datetime]
% --> wavelength - Corresponding wavelength of data, in [nano meters] 
% --> file       - file name [string]
% --> url        - url of the files [string]
% --> year       - year [double]
% --> month      - month [double]
% --> day        - day [double]
% --> date       - date [yyyymmdd]  in string

% Created on: 13th Aug 2019
% Created by: Nithin Sivadas, Boston University

if nargin <2
    timeMaxStr = '20 Jan 2018';
end

if nargin <1
    timeMinStr = '31 Dec 2017';
end

if datenum(timeMinStr)>datenum(timeMaxStr)
    error('Check inputs: the end time is smaller than the start time.');
end

host = 'optics.gi.alaska.edu';
remoteStoreLink = '/PKR/DASC/RAW/';
date1 = datetime(datenum(timeMinStr),'ConvertFrom','datenum');
date2 = datetime(datenum(timeMaxStr),'ConvertFrom','datenum');

dasc=ftp(host);
ME = [];

data.time = [];
data.wavelength = [];
data.file = [];

yearArr = year(date1):year(date2);

for iYear = 1:1:length(yearArr)
    try
        % Accessing folders within a particular year
        cd(dasc,[remoteStoreLink,num2str(yearArr(iYear))]);
        remoteFileList = dir(dasc);
        remoteFileListName = deblank(string(char(remoteFileList.name)));
        
        dayArr = datetime(remoteFileListName,'InputFormat','uuuuMMdd');
        date11 = datetime(year(date1),month(date1),day(date1));
        date21 = datetime(year(date2),month(date2),day(date2));
        
        % Flagging the folders that are within the time range
        indxWantedDays = dayArr>=date11 & dayArr<=date21;
        indx = 1:1:length(dayArr);
        indxNum = indx(indxWantedDays);
        
        % Going into the folders of interest
        for i = indxNum
           
           try
               
           % Collecting the file names within the folder    
           tempFileList = dir(dasc,...
                deblank(strcat(remoteStoreLink,num2str(yearArr(iYear)),'/',...
                remoteFileListName(i))));
            tempFileListName = string(char(tempFileList.name));
            
            % Extracting the time stamp and wavelength
            [timeStamp,wavelength] = fitsfiletimestamp(tempFileListName);
            timeStamp = datetime(timeStamp,'ConvertFrom','datenum');
            
            % Flagging the files that are within the time range
            indxWantedTimes = timeStamp>=date1 & timeStamp<=date2;
            
            % Collecting meta data into a structure
            data.file = [data.file ; tempFileListName(indxWantedTimes)];
            data.time = [data.time ; timeStamp(indxWantedTimes)];
            data.wavelength = [data.wavelength ; wavelength(indxWantedTimes)'];
            catch ME
            end
        end
       
    catch ME
          
    end
      
end        

close(dasc);    

% Calculating more data from time-stamp and wavelength
data.date = string(datestr(data.time,'yyyymmdd'));
data.year = year(data.time);
data.month = month(data.time);
data.day = day(data.time);

data.url = strcat('ftp://',host,remoteStoreLink,...
    num2str(data.year),'/',data.date,'/',data.file);

% If you need to debug, the last error is stored here. 
error = ME;

end

