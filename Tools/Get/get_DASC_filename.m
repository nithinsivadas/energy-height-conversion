function [fileStrName] = get_DASC_filename(dirName,timeStr)
%get_DASC_filename.m Identifies the file closest to the specified time
%in 'timeStr'
%   Detailed explanation goes here

%% Identifying all files in folder
fileStr = get_files_in_folder(dirName);

%% Generating time stamp of all files
aldtnum = fitsfiletimestamp(fileStr); %from John Swoboda
timeASI = (aldtnum-datenum('jan-01-1970'))*(24*3600);
timeASI = unix_to_matlab_time(timeASI);

%% Identifying file
timeIndex = find_time(timeASI,timeStr);
fileStrName = strcat(dirName,'\',(fileStr(timeIndex)));
end

