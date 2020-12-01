function [data, err] = get_DASC_FITS_times(timeMinStr, timeMaxStr)
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
% Updated on: 1 Dec 2020, used the curl command instead of ftp - to avoid
% connection issues. 
% Created by: Nithin Sivadas, Boston University


if nargin <2
    timeMaxStr = '20 Jan 2018';
end

if nargin <1
    timeMinStr = '31 Dec 2017';
end

if datenum(timeMinStr)>datenum(timeMaxStr)
    err('Check inputs: the end time is smaller than the start time.');
end

host = 'ftp://optics.gi.alaska.edu';
curlcmd = 'curl -ls '; %Silently lists the contents of the directory that is passed
remoteStoreLink = '/PKR/DASC/RAW/';
date1 = datetime(datenum(timeMinStr),'ConvertFrom','datenum');
date2 = datetime(datenum(timeMaxStr),'ConvertFrom','datenum');

ME = [];

data.time = [];
data.wavelength = [];
data.file = [];
err.subdirectory.name =[];
err.subdirectory.message =[];
err.directory.name =[];
err.directory.message =[];

yearArr = year(date1):year(date2);

for iYear = 1:1:length(yearArr)
    try
        % Accessing folders within a particular year
        [status,out] = system([curlcmd,host,remoteStoreLink,num2str(yearArr(iYear)),'/']);
        
        if status==0 && ~isempty(out)
        remoteFileListName = deblank(string(strsplit(out(1:end-1))))';
        
        dayArr = datetime(remoteFileListName,'InputFormat','uuuuMMdd');
        date11 = datetime(year(date1),month(date1),day(date1));
        date21 = datetime(year(date2),month(date2),day(date2));
        
        % Flagging the folders that are within the time range
        indxWantedDays = dayArr>=date11 & dayArr<=date21;
        indx = 1:1:length(dayArr);
        indxNum = indx(indxWantedDays);
        
        % Going into the folders of interest
            for i = indxNum

               % Collecting the file names within the folder 
               [statusFile,out] = system(strjoin([curlcmd,host,...
                   deblank(strcat(remoteStoreLink,num2str(yearArr(iYear)),...
                   '/',remoteFileListName(i))),'/'],''));

                    if statusFile==0 && ~isempty(out)
                    tempFileListName = string(strsplit(out(1:end-1)))';

                    % Extracting the time stamp and wavelength
                    [timeStamp,wavelength] = fitsfiletimestamp(tempFileListName);
                    timeStamp = datetime(timeStamp,'ConvertFrom','datenum');

                    % Flagging the files that are within the time range
                    indxWantedTimes = timeStamp>=date1 & timeStamp<=date2;

                    % Collecting meta data into a structure
                    data.file = [data.file ; tempFileListName(indxWantedTimes)];
                    data.time = [data.time ; timeStamp(indxWantedTimes)];
                    data.wavelength = [data.wavelength ; wavelength(indxWantedTimes)'];
                    else
                    err.subdirectory.name = [err.subdirectory.name;...
                        [curlcmd,host,...
                        deblank(strcat(remoteStoreLink,num2str(yearArr(iYear)),'/',...
                        remoteFileListName(i)))]];
                    err.subdirectory.message = [err.subdirectory.message; out];
                    end

            end
        
        else
            err.directory.name = [err.directory.name;...
                [curlcmd,host,remoteStoreLink,num2str(yearArr(iYear))]];
            err.directory.message = [err.directory.message; out];
        end
       
    catch ME
  
    end
      
end        


% Calculating more data from time-stamp and wavelength
data.date = string(datestr(data.time,'yyyymmdd'));
data.year = year(data.time);
data.month = month(data.time);
data.day = day(data.time);

data.url = strcat(host,remoteStoreLink,...
    num2str(data.year),'/',data.date,'/',data.file);

% If you need to debug, the last err is stored here. 
err.program.message = ME.message;

end

