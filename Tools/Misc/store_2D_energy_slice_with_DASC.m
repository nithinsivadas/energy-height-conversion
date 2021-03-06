function store_2D_energy_slice_with_DASC...
    (paths,timeLimStr,DASCcalFile,latLim,lonLim,amisrData,dataInv,magcoords)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Input:
% paths.fitsDirName -> contains the folder with fits files
% paths.imageDirName -> the folder where the figures will be stored
% timeLimStr.min -> '01 Mar 2011 9:00:00'
% timeLimStr.max -> '01 Mar 2011 10:00:00'
% latLim -> [63 67]
% DASCcalFile.az - azimuth calibration file of DASC
% DASCcalFile.el - elevation calibration file of DASC

fileStr = get_files_in_folder(paths.fitsDirName);

%% Generating time stamp for all FIT files
aldtnum = fitsfiletimestamp(fileStr);
timeASI = (aldtnum-datenum('jan-01-1970'))*(24*3600);
timeASI = unix_to_matlab_time(timeASI);

itimeStart = find_time(timeASI,timeLimStr.min);
itimeEnd = find_time(timeASI,timeLimStr.max);
mkdir(paths.fitsDirName,paths.imageDirName);
set(0,'DefaultFigureVisible','off'); 
hWait = waitbar(0);
for timeIndex = itimeStart:1:itimeEnd
    custom_waitbar(hWait,timeIndex-itimeStart+1,itimeEnd-timeIndex+1,...
        'Storing Images');
    figureHandle = figure; 
    plot_2D_energy_slice_with_DASC(figureHandle,datestr(timeASI(timeIndex)),...
    dataInv,amisrData,magcoords,DASCcalFile,...
    paths.fitsDirName,latLim,lonLim,paths.imageDirName);

end
delete(hWait);
set(0,'DefaultFigureVisible','on');
end

