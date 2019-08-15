function [data, error] = get_DASC_FITS_times(timeMinStr, ...
    timeMaxStr)
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
% T = table('Year','Date','File','Month','Day','Wavelength','Time');
for iYear = 1:1:length(yearArr)
    try
        cd(dasc,[remoteStoreLink,num2str(yearArr(iYear))]);
        remoteFileList = dir(dasc);
        remoteFileListName = string(char(remoteFileList.name));
        
        dayArr = datetime(remoteFileListName,'InputFormat','uuuuMMdd');
        date11 = datetime(year(date1),month(date1),day(date1));
        date21 = datetime(year(date2),month(date2),day(date2));
        indxWantedDays = dayArr>=date11 & dayArr<=date21;
        indx = 1:1:length(dayArr);
        indxNum = indx(indxWantedDays);
        for i = indxNum
           try
           tempFileList = dir(dasc,...
                strcat(remoteStoreLink,num2str(yearArr(iYear)),'/',...
                remoteFileListName(i)));
            tempFileListName = string(char(tempFileList.name));
            [timeStamp,wavelength] = fitsfiletimestamp(tempFileListName);
            timeStamp = datetime(timeStamp,'ConvertFrom','datenum');
            indxWantedTimes = timeStamp>=date1 & timeStamp<=date2;
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

data.date = string(datestr(data.time,'yyyymmdd'));
data.year = year(data.time);
data.month = month(data.time);
data.day = day(data.time);
data.url = strcat('ftp://',host,remoteStoreLink,...
    num2str(data.year),'/',data.date,'/',data.file);

error = ME;

end

