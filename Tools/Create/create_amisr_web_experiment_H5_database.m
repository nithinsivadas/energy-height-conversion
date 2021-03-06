function [status] = create_amisr_web_experiment_H5_database...
    (timeMinStr,timeMaxStr,instrumentID,outputH5File,writeModeStr)
%create_amisr_web_experiment_H5_database.h5 Downloads from amisr.com the data
% regarding amisr experiments, and stores it in an HDF5 file
%-------------------------------------------------------------------------
% Input
%------
%       timeMinStr   - Start time string
%       timeMaxStr   - Stop time string
%       instrumentID - The number of the instrument. It is 61 for Poker Flat ISR
%       outputH5File - Output database file name
%       writeModeStr -  NA
%----------------------------------------------------------------------------
% Output
%---------
%       status - 'Success' or 'Failed'
%----------------------------------------------------------------------------
% Modified: 30th May 2019
% Created : Unknown
% Author  : Nithin Sivadas
% Ref     :
% Comments:
%
%----------------------------------------------------------------------------

if nargin<4
    outputH5File = 'amisrWebDatabase.h5';
end

if  nargin<3 || isempty(instrumentID)
    instrumentID=61;
end
if nargin<2 || isempty(timeMaxStr)
    timeMaxStr = '5 Feb 2017';
end
if nargin <1 || isempty(timeMinStr)
    timeMinStr = '5 Jan 2017';
end

[table]=get_amisr_experiment_times...
(timeMinStr, timeMaxStr, instrumentID);

if instrumentID==61
    instrumentStr = 'PFISR';
end
status = 'failed';

table(:,3)=num2cell(posixtime(datetime(datenum(table(:,3)),'ConvertFrom','datenum')));
table(:,4)=num2cell(posixtime(datetime(datenum(table(:,4)),'ConvertFrom','datenum')));

try h5create(outputH5File,['/',instrumentStr,'/expId'],size(table(:,1)'));catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/expName'],size(table(:,2)'));catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/startTime'],size(table(:,3)'));catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/endTime'],size(table(:,4)')); catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/processingStatus'],size(table(:,5)')); ME; catch end

hdf5write(outputH5File,['/',instrumentStr,'/expId'],((table(:,1))'),...
['/',instrumentStr,'/expName'],(table(:,2))',...
['/',instrumentStr,'/startTime'],(table(:,3))',...
['/',instrumentStr,'/endTime'],(table(:,4))',...
['/',instrumentStr,'/processingStatus'],(table(:,5))');

status = 'Success';

end
