function [status, error] = create_dasc_H5_database...
    (timeMinStr,timeMaxStr,outputH5File,writeModeStr)
%create_amisr_web_experiment_H5_database.h5 Downloads from DASC FTP the data
% regarding DASC experiments, and stores it in an HDF5 file
%-------------------------------------------------------------------------
% Input
%------
%       timeMinStr   - Start time string
%       timeMaxStr   - Stop time string
%       outputH5File - Output database file name
%       writeModeStr -  NA
%----------------------------------------------------------------------------
% Output
%---------
%       status - 'Success' or 'Failed'
%----------------------------------------------------------------------------
% Modified: 15th Aug 2019
% Created : Unknown
% Author  : Nithin Sivadas
% Ref     :
% Comments:
%
%----------------------------------------------------------------------------

if nargin<3
    outputH5File = 'dascDatabase.h5';
end

if nargin<2 || isempty(timeMaxStr)
    timeMaxStr = '20 Jan 2018';
end
if nargin <1 || isempty(timeMinStr)
    timeMinStr = '31 Dec 2017';
end

status = 'failed';

[data, error] = get_DASC_FITS_times(timeMinStr, timeMaxStr);

% Can be edited to include more cameras
instrumentStr = 'DASC';

data.time=posixtime(data.time);

try h5create(outputH5File,['/',instrumentStr,'/time'],size(data.time'));catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/wavelength'],size(data.wavelength'));catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/file'],size(data.file'));catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/date'],size(data.date'));catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/url'],size(data.url'));catch ME; end

hdf5write(outputH5File,['/',instrumentStr,'/time'],(data.time'),...
['/',instrumentStr,'/wavelength'],(data.wavelength)',...
['/',instrumentStr,'/file'],(data.file)',...
['/',instrumentStr,'/date'],(data.date)',...
['/',instrumentStr,'/url'],(data.url)');

status = 'Success';

end
