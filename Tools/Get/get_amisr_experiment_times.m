function [table, tableHeader] = get_amisr_experiment_times...
    (timeMinStr, timeMaxStr, instrumentID)
%% get_amisr_experiment_times.m Scrapes experiment details from amisr website
%--------------------------------------------------------------------------
% Input
%------
% timeMinStr      : [String] The start time from which you want experiment
%                   details
%                   Default:'5 Jan 2017'
% timeMaxStr      : [String] The start time from which you want experiment
%                   details
%                   Default:'5 Feb 2017'
% InstrumentID    : [Integer] Madrigal ID number of the instrument
%                   Default:61
%                   This specified poker flat
% --------------------------------------------------------------------------
% Output
%-------
% table            : Cell array with details of the experiments from amisr
%                    webpage
% tableHeader      : Cell array with column names
%--------------------------------------------------------------------------
% Modified: 26th Feb 2017 
% Created : 23rd Feb 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------
%% Initializing
if  nargin<3
    instrumentID=61;
end
if nargin<2
    timeMaxStr = '5 Feb 2017';
end
if nargin <1
    timeMinStr = '5 Jan 2017';
end

if datenum(timeMinStr)>datenum(timeMaxStr)
    error('Check inputs: the end time is smaller than the start time.');
end

%% Populating Table
tables.idTableBy.plaintextPreceedingTable='Filter'; 
nMonths = months(timeMinStr,timeMaxStr) +1;
startMonth = str2num(datestr(timeMinStr,'mm'));
startYear = str2num(datestr(timeMinStr,'yyyy'));
table = [];

for iMonth = 1:1:nMonths
    try 
    mm = num2str(startMonth);
    yyyy = num2str(startYear);
    % Link to the webpage with experiment information
    amisrDataLink = ['http://amisr.com/database/'...
        ,num2str(instrumentID),'/exp-list/',yyyy,'/',...
        mm,'/'];
    tempTableCell = htmlTableToCell(amisrDataLink,tables);
    if iMonth==1
        tableHeader = tempTableCell(1,:);
        timeIndxBeg=find_time(datenum(tempTableCell(2:end,3)),timeMinStr)+1;
    else 
        timeIndx = 2;
    end
    if iMonth==nMonths
        % Finding indx that includes the timeMaxStr day
        timeIndxEnd=find_time(datenum(tempTableCell(2:end,3))-1,timeMaxStr);
    else
        timeIndxEnd=size(tempTableCell,1);
    end
    % Appending cell array with data
    table = [table ; tempTableCell(timeIndxBeg:timeIndxEnd,:)];
    catch ME
    end
    if startMonth<12
        startMonth = startMonth+1;
    else 
        startMonth = 1;
        startYear = startYear + 1;
    end
    
end


end

