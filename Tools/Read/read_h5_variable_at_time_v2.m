function [varValue] = read_h5_variable_at_time_v2(h5FileStr, datasetPath,...
    thisTimeIndx, thisTimeStr)
%READ_H5_VARIABLE_AT_TIME_V2 Reads variable at timeIndx or timeStr

%----------------------------------------------------------------------
    % Input
    %--------     
    %   h5FileStr   - The path of the hdf5 file
    %   datasetPath - Internal dataset path including its name
    %   thisTimeIndx- Index of the time of interest
    %   thisTimeStr - Time string in format 'dd-mm-yyyy HH:MM:SS.FFF'
    %----------------------------------------------------------------------
    % Output
    %--------
    %   varValue    - Data within the variable for that time instant
    %----------------------------------------------------------------------
    % Comments
    % Provide timeIndx for faster access
    % If both timeIndx and timeStr are provided, then timeIndx will get
    % preference
    % If datasets with no time variation are requested, the dataset will be
    % read without any change, and presented as output.
    %
    % Last Updated: 23rd Nov 2018,30th Sep 2018
    % Create by   : Nithin Sivadas


if nargin < 4
    thisTimeStr = [];
end

if nargin < 3
    error('Please specify either thisTimeIndx, or thisTimeStr');
end

temp = strsplit(datasetPath,'/');
groupStr = strjoin(temp(1:end-1),'/');
timeStr = strcat(groupStr,'/time');
[status, info, ~] = ish5dataset(h5FileStr, datasetPath);
if ~status
    error(['Dataset ',datasetPath,' incorrect or does not exist']);
end
% Getting thisTimeIndx if not provided (might take a bit longer)
if isempty(thisTimeIndx)
    if isempty(thisTimeStr)
        error('Please specify either thisTimeIndx, or thisTimeStr');
    else

        [statusTime, infoTime] = ish5dataset(h5FileStr, timeStr);
        if statusTime
            try
                units=h5readatt(h5FileStr,timeStr,'Units');
            catch
                units = '[matlab units]';
            end
            time = h5read(h5FileStr,timeStr);
        else
            warning(['No time array for the group,',...
                'using the default time array: /energyFluxFromMaxEnt/time']);
            timeStr = '/energyFluxFromMaxEnt/time';
            [statusTime, infoTime] = ish5dataset(h5FileStr, timeStr);
            if statusTime
                try
                    units=h5readatt(h5FileStr,timeStr,'Units');
                catch
                    units= '[matlab units]';
                end
                time = h5read(h5FileStr,timeStr); %This is currently matlab time :/
            else
                error('Default Time array /energyFluxFromMaxEnt/time does not exists');
            end
        end
        if ~strcmp(units,'[matlab units]')
            time = unixtime2matlab(time);
        end
        thisTimeIndx = find_time(time,thisTimeStr);
    end
    nTime = length(time);
else
    [statusTime, infoTime] = ish5dataset(h5FileStr, timeStr);
    if ~statusTime
        warning(['No time array for the group,',...
                'using the default time array: /energyFluxFromMaxEnt/time']);
        timeStr = '/energyFluxFromMaxEnt/time';
        [statusTime, infoTime] = ish5dataset(h5FileStr, timeStr);
        if ~statusTime
            error('Default Time array /energyFluxFromMaxEnt/time does not exists');
        end
    end
    nTime = max(infoTime.Dataspace.Size);
end



dataSize = info.Dataspace.Size;
timeDim = find(dataSize==nTime);
nDim = length(dataSize);

startIndx = ones(1,nDim);
countIndx = dataSize;
if ~isempty(timeDim)
    if thisTimeIndx <= countIndx(timeDim)
        countIndx(timeDim) = 1;
        startIndx(timeDim) = thisTimeIndx;
    else
        error('Time Index exceeds the size of records in h5 file');
    end
end

varValue = squeeze(permute(h5read(...
    h5FileStr,datasetPath,startIndx,countIndx...
    ),fliplr(1:1:nDim)));


end
